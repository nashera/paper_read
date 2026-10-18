# 08 癌种预测与外部 BRCA 验证

> **目标文献**：Tang Z, et al. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* Nature Communications (2026) 17:10800. DOI: 10.1038/s41467-026-77538-5
> **模块对应原文**：Fig. 3a–f；Supplementary Figs. 4–5；Results「Prediction of cancer types across primary and metastatic samples」
> **前置依赖**：`05`（多分类评估）、`07`（训练好的 Stage 2 模型）；`raw/brca_cellxgene/`
> **复现难度 / 预计工时 / 计算资源**：★★★☆☆ / 2–3 人日 / V100 推理数小时

---

## ① 科学问题与逻辑位置

- **上游**：Stage 2 多分类模型 + 泛癌测试集 + 外部 BRCA。
- **下游**：`figures/Fig3*`；13 号 QC 门槛；06 号的泛化结论。
- **本模块回答**：**"模型能否在未见过的 donor、甚至未见过的数据集上，判断肿瘤的原发部位？"**

这对 **CUP（原发灶不明癌）** 有直接临床意义：转移灶可活检，但原发灶找不到；如果表达谱能预测组织来源，就能指导治疗。

---

## ② 背景原理

### 2.1 两阶段推理

- Stage 1 判 tumor/normal；
- 对 tumor bag，Stage 2 输出 **29 个癌种标签**的概率；`2-bag_eval.py` 再把 Stage-1 判 normal 的 bag 置为 `Y_hat=29`，形成第 30 个**预测桶（normal/rejection）**（D16）。
- sample 级：对 sub-sample 的概率**投票/求和**（05 号）。

### 2.2 原发 vs 转移

- **原发（primary）**：组织来源明确，表达谱保留起源特征。
- **转移（metastatic）**：微环境改变、上皮-间质转化（EMT）、器官特异基质信号，可能"冲淡"起源特征。
- 原文的 Sankey 图分别展示 primary/metastatic 的对角线匹配度。

### 2.3 外部 BRCA 的"单类"评估

BRCA 数据集全部是乳腺癌，**没有负类/其他癌种**，所以：

- 不能算多分类 AUC/混淆矩阵；
- 能算的是"被判为 BRCA 的比例"（accuracy-like）；
- 以及"被误判成哪些癌种"（误分类分布）。
- 原文报告"almost all primary samples correctly identified as BRCA, only a small fraction of metastatic samples misclassified"。

### 2.4 误分类的生物学解释（原文 Supplementary Fig. 5）

| 误判 | 可能原因 |
|---|---|
| → THCA | 高淋巴浸润造成"信号稀释"（thyroid 与 lymphoid 的表达重叠） |
| → SARC-EWING | 骨转移的基质信号（Ewing 肉瘤常位于骨） |
| TNBC → PAAD/OV | 三阴性乳腺癌与高级别胰腺/卵巢癌的分子相似性 |

> ⚠️ **导师批注**：这些是**事后解释**，不是因果。复现时应把它们当作"假设生成"，并明确写"原文推测"。不要写成"证明了三阴性乳腺癌与卵巢癌同源"。

### 2.5 Bag embedding 的 UMAP

- 取 Stage 1 的 bag embedding（256 维）；
- 用 scanpy 的 neighbors + UMAP 可视化；
- 颜色：癌种（生物保守）、数据集来源（批次）；
- 期望：癌种分簇、数据集混合（scIB 量化）。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "scTransMIL achieved better classification accuracy for both primary and metastatic tumors" | 原发/转移均好 | 多分类指标 | 未给数值 | 以 Source Data 为准 |
| "the model's robust performance on metastatic specimens is especially noteworthy" | 转移泛化 | 分层报告 | — | — |
| "BRCA dataset … 400,000 cells from 37 samples … all bags assigned BRCA" | 外部单类 | 推理 | 给了规模 | — |
| "almost all primary samples correctly identified as BRCA, and only a small fraction of metastatic samples misclassified" | 误分类少 | 比例 | 未给具体 | 以图为准 |
| "Thyroid Carcinoma (THCA) misclassifications appear to be linked to 'signal dilution' from high lymphoid infiltration" | 误判解释 | 事后分析 | — | 标为推测 |
| "samples with a triple-negative status were frequently misclassified as PAAD or OV" | TNBC 误判 | 事后分析 | — | 标为推测 |

---

## ④ 逐步可运行代码

### 4.1 多分类评估

```python
# code/08_01_multicls_eval.py
import pandas as pd, numpy as np
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix

df = pd.read_csv("results/metrics_multicls_scTransMIL_fold0.csv")
acc = accuracy_score(df.Y, df.Y_hat)
f1 = f1_score(df.Y, df.Y_hat, average="macro", zero_division=1)
cm = confusion_matrix(df.Y, df.Y_hat)
print("acc", acc, "macro-F1", f1)
pd.DataFrame(cm).to_csv("results/08_confusion_multicls.csv", index=False)
```

### 4.2 Sankey 数据

```python
# code/08_02_sankey.py
# 输出：source(真实癌种) -> target(预测癌种) -> value
import pandas as pd
df = pd.read_csv("results/metrics_multicls_scTransMIL_fold0.csv")
links = df.groupby(["Y", "Y_hat"]).size().reset_index(name="value")
links.to_csv("results/08_sankey_links.csv", index=False)
```

> 画 Sankey 可用 `plotly.graph_objects.Sankey`；节点 = 癌种，连线宽度 = 样本数。原发/转移分别画（需要 `is_metastatic` 列；若数据没有，从 `Tissue`/样本注释推断并在 D 表记录）。

### 4.3 BRCA 推理

```python
# code/08_03_brca_infer.py
import glob, torch, pandas as pd
from model_scTransMIL import ScTransMILCode
from torch.utils.data import DataLoader
# 复用 05 号的 BagDataset / predict_bags / aggregate_multiclass
# 1) BRCA h5ad → 02 号 gene map/tokenize（输出 processed/tokenized/brca/*.pt）
# 2) 载入 Stage-1 / Stage-2 checkpoint
stage1 = ScTransMILCode("raw/pretrained/gf-6L-30M-i2048", n_classes=2)
stage1.load_state_dict(torch.load("results/best_bicls.pt", map_location="cpu"))
stage2 = ScTransMILCode("raw/pretrained/gf-6L-30M-i2048", n_classes=29)
stage2.load_state_dict(torch.load("results/best_multicls.pt", map_location="cpu"))
# 3) 对每个 bag 预测；仅对 Stage-1 阳性的 bag 跑 Stage-2
#    （DataLoader 用 05 号 BagDataset，返回 ids/mask/label/slide_id）
# 4) 聚合到 sample；统计 predicted == BRCA 的比例
```

```python
# 关键统计
df = pd.read_csv("results/brca_predictions.csv")
print("samples:", df.sample_id.nunique())
print("predicted BRCA:", (df.pred == "BRCA").mean())
print(df[df.pred != "BRCA"].pred.value_counts())
```

### 4.4 误分类的细胞组成分析

```python
# code/08_04_misclass_composition.py
import scanpy as sc, pandas as pd
# 对每个误分类 sample，计算细胞类型比例（如 lymphoid/stromal/epithelial）
# 与正确分类 sample 比较
```

### 4.5 bag embedding UMAP

```python
# code/08_05_bag_umap.py
import scanpy as sc, numpy as np, pandas as pd
emb = np.load("processed/embeddings/bag_embedding_pancancer_brca.npy")
meta = pd.read_csv("processed/embeddings/bag_meta.csv")
ad = sc.AnnData(emb)
ad.obs = meta.set_index(ad.obs_names)
sc.pp.neighbors(ad, n_neighbors=15, use_rep="X")
sc.tl.umap(ad, random_state=20261018)
ad.write("processed/embeddings/bag_umap.h5ad")
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| 多分类类别 | 29 个癌种 + 1 个 normal 桶（D16） | — | 报告时明确桶含义 |
| 聚合 | 概率和 argmax | 投票 | 结果略变 |
| BRCA 判据 | 预测为 BRCA | top-1 | — |
| UMAP `n_neighbors` | 15 | 5–50 | 视觉分簇 |
| 误分类阈值 | 无 | — | 描述性分析 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| 原发癌种准确率 | 高（对角线） | 以 Source Data 为准 | ⚠️ |
| 转移癌种准确率 | 略低但可辨 | — | ⚠️ |
| BRCA 正确率 | 接近 1（原发） | — | ⚠️ |
| BRCA 误分类 | 少量转移 | — | ⚠️ |
| 误分类模式 | THCA/SARC-EWING/PAAD/OV | 方向一致 | ⚠️ |
| bag UMAP | 癌种分簇、数据集混合 | 定性 | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| BRCA 全判成某一类 | 类别先验/标签映射错 | 检查 Stage 2 的 label_dict |
| 多分类准确率虚高 | 测试集类别不平衡 | 报 macro-F1 与混淆矩阵 |
| Sankey 行列不对 | Y/Y_hat 编码未映射回名称 | 用 label_dict 反查 |
| 转移样本全错 | 训练集转移样本少 | 分层报告；不要合并 |
| UMAP 按数据集分 | 批次效应 | 这是 scIB 要量化的；与原文对照 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 3a–b | `figures/Fig3AB_sankey.*` | 原发/转移 |
| Fig. 3c | `figures/Fig3C_brca_sankey.*` | BRCA |
| Fig. 3d | `figures/Fig3D_accuracy.*` | 方法对比 |
| Fig. 3e–f | `figures/Fig3EF_brca_umap.*` | 细胞类型/tumor score |
| Fig. 3g–h | `figures/Fig3GH_bag_umap.*` | **由 11 号生成**；本模块只消费 |
| Fig. S4–S5 | `results/08_*` | 癌种清单/误分类 |

> **交付物**：`results/08_confusion_multicls.csv`、`results/08_sankey_links.csv`、`results/brca_predictions.csv`、`figures/Fig3*`。完成后进入 `09_单细胞可解释性与基因注意力.md`。
