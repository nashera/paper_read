# 03 模型架构：rank 编码、Transformer、门控与聚合

> **目标文献**：Tang Z, et al. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* Nature Communications (2026) 17:10800. DOI: 10.1038/s41467-026-77538-5
> **模块对应原文**：Fig. 1c–d；Methods「Instance encoder」「Gated encoder」「Bag aggregator」；作者代码 `models/model_clam.py`、`models/model_mil.py`
> **前置依赖**：`02_数据预处理与MIL样本构造.md`（`input_ids`/`attention_mask`）；`raw/pretrained/gf-6L-30M-i2048/`
> **复现难度 / 预计工时 / 计算资源**：★★★★☆ / 3–5 人日 / GPU V100 32 GB（训练时 batch=1 bag，≤256×1024 tokens）

---

## ① 科学问题与逻辑位置

- **上游**：每个 bag 的 `input_ids`（N cells × 1024 tokens）与 `attention_mask`。
- **下游**：`processed/embeddings/` 的细胞/ bag 嵌入；`results/cell_scores` 的 tumor score；`results/attention_*` 的基因注意力。
- **本模块回答**：**"一个没有细胞级标签的模型，怎样从 token 序列算出细胞嵌入、细胞级分数和 bag 级诊断？"**

---

## ② 背景原理

### 2.1 形式化

设 bag $j$ 有 $m_j$ 个 instance（细胞），第 $i$ 个细胞的 token 序列为 $\mathbf t_i\in\mathbb N^{L}$（$L=1024$）。

$$
\mathbf h_i^{(0)}=\text{Transformer}(\mathbf t_i)\in\mathbb R^{d},\qquad d=256,
$$

$$
\mathbf h_i^{\text{gate}}=\tanh(\mathbf W_1\mathbf h_i^{(0)})\odot\sigma(\mathbf W_2\mathbf h_i^{(0)}),
$$

$$
a_i=\mathbf w^\top\mathbf h_i^{\text{gate}}\in\mathbb R,
$$

$$
\mathbf h_j^{\text{bag}}=\sum_{i=1}^{m_j}a_i\,\mathbf h_i^{\text{gate}},
$$

$$
\hat y_j=\operatorname{softmax}(\mathbf W_c\mathbf h_j^{\text{bag}}).
$$

> 这是**正文 Eq.1–3** 的写法。作者代码的实现略有不同（§2.5），必须两套都理解。

### 2.2 自注意力为什么适合 rank 序列

Transformer 的一层：

$$
\mathbf Q=\mathbf X\mathbf W_Q,\ \mathbf K=\mathbf X\mathbf W_K,\ \mathbf V=\mathbf X\mathbf W_V,
$$

$$
\operatorname{Attn}(\mathbf Q,\mathbf K,\mathbf V)=\operatorname{softmax}\!\left(\frac{\mathbf Q\mathbf K^\top}{\sqrt{d_k}}\right)\mathbf V.
$$

- rank 序列里，位置 = 表达名次；注意力可以学"前 50 名基因与后 500 名基因的交互"。
- 多头注意力允许不同头关注不同"基因模块"（如增殖、免疫、分化）。
- 与 MLP 不同，Transformer 对**输入顺序**敏感（需要位置编码）——Geneformer 使用 rank 位置编码，名次本身携带信息。

### 2.3 实例编码器的两种实现

| 实现 | 论文描述 | 代码 |
|---|---|---|
| 池化 | "averaging the embeddings of each gene … 256-dim cell embedding" | `BertModel(...).pooler_output`（CLS token 经 dense+tanh） |
| 层数 | 6 层 Transformer | Geneformer `gf-6L-30M-i2048`（6 层） |
| 维度 | 256 | 以 `config.json` 为准；作者 `CLAM_SB` 的 `L=256` 说明 pooler 输出 256 |

> ⚠️ **导师批注**：这是 D3 的第一处"正文 vs 代码"差异。均值池化与 CLS 池化在 rank 序列上的行为不同：均值池化更"平滑"，CLS 池化更"任务驱动"。**以代码为准复现**，另写一版 mean pooling 做对照。

### 2.4 门控编码器

正文：

$$
\mathbf h_1=\tanh(\mathbf W_1\mathbf h_{\text{init}}),\qquad
\mathbf h_2=\sigma(\mathbf W_2\mathbf h_{\text{init}}),\qquad
\mathbf h_{\text{gate}}=\mathbf h_1\odot\mathbf h_2,
$$

$$
a=\operatorname{linear}(\mathbf h_{\text{gate}}).
$$

- $\mathbf h_2$ 是"软掩码"：sigmoid 输出 0–1，逐维抑制与 bag 标签无关的特征。
- $\mathbf h_1$ 是"候选特征"：tanh 输出 −1–1。
- 二者逐元素相乘 = 门控。

**代码**（`Attn_Net_Gated`）：

```python
a = Linear(L, D) + Tanh
b = Linear(L, D) + Sigmoid
A = a * b
A = Linear(D, 1)(A)      # 每个 instance 一个标量
A = softmax(A.T, dim=1)  # 在 bag 内归一化
M = A @ h                # ★ 用原始 h，不是门控后的 A_raw
```

- $L=256$（输入嵌入维），$D=64$（注意力隐藏维）。
- **注意**：代码的 `Attn_Net_Gated.forward` 返回 `(A, x)`，即聚合用的是**原始细胞嵌入** $h$，门控只影响注意力权重；这与正文"聚合 $\mathbf h_{\text{gate}}$"不同（D3 第二处）。
- **三个量必须分开**（D22）：
  1. `a_i`（正文 Eq.2）：门控分支的**未归一化**标量，用于 adaptive loss 排序；
  2. `attention_weight = softmax(a)`：bag 内归一化权重，用于聚合 `M = A @ h`；
  3. `cell_score`：作者 `eval_utils.py` 实际保存的是 **softmax 后**的 attention weight（用于单细胞 AUC/UMAP）。与 `a_i` 不是同一个量。

> **D3 第四处差异**：正文 Eq.6 的 instance 分类器输入是 $\mathbf h_{\text{gate}}$；作者代码 `inst_eval` 的 `all_instances` 来自**原始** $\mathbf h$（`model_clam.py:143`）。本项目以代码为准，正文公式版作为对照。

### 2.5 Bag 聚合

正文 Eq.3：$\mathbf h_{\text{bag}}=\sum_i a_i\mathbf h_i^{\text{gate}}$。

代码：$\boldsymbol\alpha=\operatorname{softmax}(\mathbf A)$，$\mathbf h_{\text{bag}}=\sum_i\alpha_i\mathbf h_i$。

**为什么代码用 softmax**：softmax 使权重和为 1、梯度稳定；未归一化的 $a_i$ 会随 bag 大小缩放，导致 bag embedding 范数不稳定。正文的写法在数学上更"像 OR 逻辑"，但训练上更抖。

### 2.6 二阶段分类头

```python
# Stage 1
logits1 = Linear(256, 2)(h_bag); y1 = softmax(logits1)
# Stage 2（仅 tumor bag）
logits2 = Linear(256, n_classes)(h_bag); y2 = softmax(logits2)
```

- **正文**：Stage 2 复用 Stage 1 的 bag embedding（"using the bag embedding from Stage 1"）。
- **代码**：`main.py` 的 `task_1_tumor_vs_normal` 与 `task_2_tumor_subtyping` 是**两个独立训练任务**，各自一个 `CLAM_SB`；`2-bag_eval.py` 只在推理时把 Stage-1 判 normal 的 bag 置为第 30 个预测桶，**没有传递 bag embedding**。以代码为准（D3）。
- 代码里二分类与多分类是**两个独立训练任务**（`task_1` / `task_2`），不是同一次前向的两阶段；正文 Fig. 1c 的"两阶段"是推理流程的串联（先用 Stage 1 判 tumor，再对 tumor bag 跑 Stage 2 模型）。D3 第三处。

### 2.7 参数量估算

| 组件 | 参数量（约） |
|---|---|
| Geneformer V1（6 层，hidden 256） | ~10 M（"30M" 指 Genecorpus-30M 细胞数，不是参数量） |
| 门控注意力（2×256→64，64→1） | ~33 K |
| 分类头（256→2 / 256→29） | ~0.5 K / ~7.5 K |
| 实例分类器（256→2，每类一个） | 与类数线性 |
| **合计** | **~10 M**（以 `config.json` 实测为准） |

> ~10 M 参数在 V100 32 GB 上、batch=1 bag（≤256×1024 tokens）可训练；显存主要被注意力矩阵占据——注意 4 个 head × 6 层的注意力分数是 $256\times4\times1024^2\times4$ B ≈ 4.3 GB/层，eager 模式下 6 层约 26 GB，需 SDPA/FlashAttention 或 gradient checkpointing 才能压到 10 GB 以下（H16）。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应实现 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "the transcriptomic profile of each cell is encoded into a ranked gene sequence" | rank 序列 | 02 号 | 未给长度 | D4 |
| "The ranked gene sequence … is fed into a 6-layer Transformer to generate initial cell embeddings." | 6 层 Transformer | Geneformer 6L | 给了层数 | D10 |
| "The initial embedding of any cell is generated by averaging the embeddings of each gene in the gene sequence, resulting in a 256-dimensional initial cell embedding." | 均值池化 + 256 维 | 代码为 CLS pooler | 给了维度 | D3 |
| "a learnable gate filters bag label-irrelevant information … generating gated embeddings, which are then linearly projected to produce instance-level tumor scores." | 门控 + tumor score | `Attn_Net_Gated` | 未给隐藏维 | D3（D=64） |
| "all gated embeddings and tumor scores in the bag are weighted and summed into a bag embedding" | 加权求和 | 代码 `A @ h` | 给了公式 | D3 |
| "If classified as a tumor bag, it enters Stage 2, which employs a multi-classifier to predict the specific cancer type using the bag embedding from Stage 1." | 两阶段 | 两个独立任务 | — | D3 |

---

## ④ 逐步可运行代码

### 4.1 按正文公式实现（教学版）

```python
# code/03_00_model_paper.py
import torch, torch.nn as nn
from transformers import BertModel

class ScTransMILPaper(nn.Module):
    """正文 Eq.1-3 的实现：mean pooling + gated embedding 聚合。"""
    def __init__(self, bert_dir, n_classes=2, d=256, gate_dim=64):
        super().__init__()
        self.encoder = BertModel.from_pretrained(bert_dir)
        self.gate1 = nn.Linear(d, gate_dim)
        self.gate2 = nn.Linear(d, gate_dim)
        self.score = nn.Linear(gate_dim, 1)
        self.cls = nn.Linear(d, n_classes)

    def forward(self, input_ids, attention_mask):
        # (N, L) -> (N, d)
        out = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        h = out.last_hidden_state.mean(dim=1)          # 正文：均值池化
        h1 = torch.tanh(self.gate1(h))
        h2 = torch.sigmoid(self.gate2(h))
        h_gate = h1 * h2                               # Eq.1
        a = self.score(h_gate)                         # Eq.2, (N,1)
        h_bag = (a * h_gate).sum(dim=0, keepdim=True)  # Eq.3
        logits = self.cls(h_bag)
        return logits, a.squeeze(1), h, h_bag
```

### 4.2 按作者代码实现（复现主路径）

```python
# code/03_01_model_code.py
import torch, torch.nn as nn, torch.nn.functional as F
from transformers import BertModel

class AttnNetGated(nn.Module):
    def __init__(self, L=256, D=64, dropout=0.0, n_classes=1):
        super().__init__()
        self.a = nn.Sequential(nn.Linear(L, D), nn.Tanh())
        self.b = nn.Sequential(nn.Linear(L, D), nn.Sigmoid())
        if dropout > 0:
            self.a.append(nn.Dropout(dropout)); self.b.append(nn.Dropout(dropout))
        self.c = nn.Linear(D, n_classes)

    def forward(self, x):
        A = self.c(self.a(x) * self.b(x))   # (N,1)
        return A, x                          # ★ 返回原始 x（D3）

class ScTransMILCode(nn.Module):
    """作者 models/model_clam.py 的 CLAM_SB 等价实现。"""
    def __init__(self, bert_dir, n_classes=2, embed_dim=256, gate_dim=64):
        super().__init__()
        self.encoder = BertModel.from_pretrained(bert_dir)
        self.attention = AttnNetGated(L=embed_dim, D=gate_dim)
        self.cls = nn.Linear(embed_dim, n_classes)

    def forward(self, input_ids, attention_mask, return_attn=False):
        h = self.encoder(input_ids=input_ids, attention_mask=attention_mask).pooler_output
        A, h = self.attention(h)             # A: (N,1)
        A = A.transpose(1, 0)                # (1,N)
        A_raw = A
        A = F.softmax(A, dim=1)              # (1,N)
        M = torch.mm(A, h)                   # (1,d) ★ 原始 h
        logits = self.cls(M)
        if return_attn:
            # 返回 per-instance h（自适应聚类损失需要），以及 bag embedding M
            return logits, A_raw, A, h, M
        return logits
```

### 4.3 与作者仓库对齐

```python
# code/03_02_use_author_model.py
# 直接调用作者 models/model_clam.py 的 CLAM_SB，确保逐层一致
import sys, torch
sys.path.insert(0, "raw/pretrained/reference_code")
from models.model_clam import CLAM_SB

model = CLAM_SB(gate=True, size_arg="small", dropout=0.25,   # 与 main.py 默认一致
                k_sample=8, n_classes=2, embed_dim=256)
# 输入：features_tuple = (input_ids, attention_mask)
logits, Y_prob, Y_hat, A_raw, results = model((ids, mask), label=torch.tensor([1]))
print(logits.shape, A_raw.shape)
```

> ⚠️ 作者 `CLAM_SB.__init__` 里 `BertModel.from_pretrained("./models/gf-6L-30M-i2048/")` 是**硬编码相对路径**；在本项目里要改成 `raw/pretrained/gf-6L-30M-i2048`，或把权重软链到 `./models/`。这是复现时第一个会报错的地方。

### 4.4 参数量与显存检查

```python
# code/03_03_model_check.py
# 依赖 model_scTransMIL.py（03_01 的内容另存）
from model_scTransMIL import ScTransMILCode
model = ScTransMILCode("raw/pretrained/gf-6L-30M-i2048", n_classes=2)
n = sum(p.numel() for p in model.parameters())
print("params:", n / 1e6, "M")

ids = torch.randint(0, 20000, (256, 1024)).cuda()
mask = torch.ones_like(ids).cuda()
with torch.no_grad():
    logits = model(ids, mask)
print(logits.shape)
print("peak mem (MB):", torch.cuda.max_memory_allocated() / 1024**2)
```

**预期**：~10 M 参数；256×1024 输入在 V100 上的前向峰值显存**取决于 attention 实现**：SDPA/FlashAttention 约 2–4 GB，eager 约 26 GB。运行时记录 `attn_implementation` 与实测值。

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| 池化 | CLS（代码）/ mean（正文） | CLS/mean | 细胞嵌入分布不同（D3） |
| 聚合 | softmax(A) @ h（代码） | Σ a h_gate（正文） | bag 嵌入尺度/注意力锐度不同 |
| `embed_dim` | 256 | 128–512 | 与 Geneformer hidden 必须一致 |
| `gate_dim` | 64 | 32–128 | 注意力容量 |
| `dropout` | 0.25（`main.py` 默认；直接实例化 `CLAM_SB` 的默认是 0.0） | 0.0–0.5 | 小数据易过拟合 |
| 分类头 | 单层线性 | MLP | 多分类容量 |
| `max_len` | 1024 | 512–2048 | 显存与信息量（D4） |
| 是否冻结 encoder | 训练时全量更新 | 冻结/分层学习率 | 小数据微调时冻结更稳（D13） |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期 | 允许波动 | Go/No-Go |
|---|---|---|---|
| 参数量 | ~10 M | 9–11 M | ⚠️ |
| 前向输出 | logits (1,2)，A (1,N) | 必须 | ❌ |
| 注意力权重和 | =1（softmax） | 1e-5 | ✅ |
| 细胞嵌入维度 | 256 | 必须 | ❌ |
| 与作者 `CLAM_SB` 输出 | 同输入下 logits 接近 | <1e-4（权重同） | ❌ |
| 显存（前向） | 实测 + 记录 attention 实现；SDPA 2–4 GB / eager ~26 GB | — | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| `from_pretrained` 找不到 | 权重路径硬编码 | 改路径或软链 |
| hidden size 不匹配 | embed_dim 与 config 不一致 | 打印 `config.hidden_size` |
| 注意力形状错 | `A.transpose` 维度 | 保持 (N,1)→(1,N) |
| 全为 NaN | softmax 前有 inf / mask 未用 | 检查 `attention_mask` 与 padding |
| 显存 OOM | 1024 tokens × 256 cells | 降 `max_len` 或 bag |
| 与作者输出不一致 | 池化/聚合实现不同 | 用 4.3 直接调作者模型对照 |
| `pooler_output` 为 None | 模型 config 无 pooler | 改用 `last_hidden_state[:,0]` |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 1c（两阶段流程） | `ScTransMILCode` 前向 | Stage 1/2 两个模型 |
| Fig. 1d（门控与聚合） | `AttnNetGated` + `A@h` | 与正文 Eq.1–3 的差异见 D3 |
| Fig. S1（架构细节） | `03_03_model_check.py` | 参数量/显存 |
| Methods「Instance/Gated encoder」「Bag aggregator」 | 本文件 §2 | — |

> **交付物**：`code/03_01_model_code.py`、`code/03_03_model_check.py`、`logs/03_model_summary.txt`。完成后进入 `04_损失函数与训练策略.md`。
