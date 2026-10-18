# 10 AML 微调、亚型与生存

> **目标文献**：Tang Z, et al. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* Nature Communications (2026) 17:10800. DOI: 10.1038/s41467-026-77538-5
> **模块对应原文**：Fig. 5a–f；Supplementary Figs. 8–12；Results「Deciphering novel AML heterogeneity and subtypes with scTransMIL」
> **前置依赖**：`04`（训练/微调）、`05`（推理）、`09`（tumor score/attention）；`raw/aml_GSE116256/`、`raw/bulk_aml_dbgap/`
> **复现难度 / 预计工时 / 计算资源**：★★★★★ / 5–8 人日 / V100 微调数小时 + bulk 生存分钟级；dbGaP 申请数周

---

## ① 科学问题与逻辑位置

- **上游**：泛癌预训练模型；AML 数据集（新病种，不在 30 类中）。
- **下游**：`results/aml_subtypes.csv`、`results/survival_aml.csv`、`figures/Fig5*`。
- **本模块回答**：**"一个在 30 类癌种上训练的模型，能否只用几十个样本级标签，适配到一个全新的疾病，并从中发现临床有意义的亚型？"**

---

## ② 背景原理

### 2.1 AML 的生物学

- AML 是髓系造血干/祖细胞的恶性克隆；
- 分化阻滞在不同阶段 → 形态/表达谱异质；
- 正常造血分化轨迹：HSC → MPP → GMP/Prog → ProMono → Mono → cDC；
- 恶性细胞常"卡"在某个分化状态（HSC-like、GMP-like、Mono-like…）；
- HSC-like 白血病干细胞（LSC）通常耐药、易复发 → 生存差。

### 2.2 微调策略

原文比较三种：

| 策略 | 冻结 | 训练 | 结果 |
|---|---|---|---|
| 不微调 | 全部 | 无 | 与 AML 不匹配 |
| 只训 aggregator | encoder | 分类头/聚合 | 次优 |
| 最后 2 个 Transformer 块 + aggregator | 其余 | 部分 encoder | **最优** |

> 正文说"selectively fine-tuning the last two transformer blocks alongside the Bag aggregator … yielded optimal performance"。代码里需确认是否实现"最后 2 块"的开关（D13）。

### 2.3 Donor setting 的微调标签

- AML 数据有 21 donors（5 healthy + 16 AML）；
- 微调只用 **20 tumor + 3 normal samples** 的 donor 级标签；
- 测试 **15 tumor + 3 normal**；
- Donor setting：同一 donor 的所有 sub-sample 用同一标签（包括可能正常的 sub-sample）→ 标签噪声。

### 2.4 标签效率

- down-sampling：10 tumor + 3 normal 时性能进入平台；
- 说明模型不需要大量标注，适合罕见病。

### 2.5 亚型发现

1. 用微调后的模型给每个细胞打分；
2. **tumor score > 0.2** 的细胞视为恶性（保守阈值）；
3. 在 **instance embedding** 上做 kNN，把恶性细胞对齐到**正常分化状态**（HSC-like/GMP-like/Prog-like/ProMono-like/Mono-like/cDC-like）；
4. 每个亚型 >1000 个肿瘤细胞；
5. 轨迹：HSC-like → GMP-like → Mono-like。

### 2.6 生存验证

- bulk AML 671 例（dbGaP）；
- 对每个亚型取 attention top-100 基因中的 **top 30**；
- 对每个患者计算这 30 个基因的平均表达；
- 每个亚型取平均表达最高的 **top 50 患者**（互不重叠）→ 6×50=300；
- 6 组质心层次聚类 → 与分化层级一致；
- KM 比较各亚型；HSC-like vs GMP-like 的 log-rank。

> ⚠️ **导师批注**：这是"用单细胞定义的亚型给 bulk 患者分组"的典型做法，但**分组阈值（top 50/top 30）是人为的**，会引入选择偏差。复现时必须做敏感性分析（top 30/50/100；top 20/50 基因），并报告"分组是否稳定"。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "AML dataset comprised 35 tumor and 6 normal samples" | 35+6 | 数据 | 给了 | — |
| "20 tumor and 3 normal samples from 21 donors for fine-tuning and the remaining 15 tumor and 3 normal samples for evaluating" | 20+3 / 15+3 | 划分 | 给了 | D13 |
| "we compared three experimental settings" | 三种微调 | 策略 | 给了 | D13 |
| "selectively fine-tuning the last two transformer blocks alongside the Bag aggregator … yielded optimal performance" | 最后 2 块 + aggregator | 冻结策略 | 给了 | D13 |
| "performance … plateaus with as few as 10 tumor and 3 normal training samples" | 10+3 平台 | 下采样 | 给了 | — |
| "conservative tumor score threshold of >0.2 followed by kNN-based annotation" | 0.2 + kNN | 亚型 | 给了 | D14 |
| "six distinct AML malignant subtypes (each comprising >1000 tumor cells)" | 6 亚型 | 聚类 | 给了 | D14 |
| "top 50 patients … top 30 attention-prioritized signature genes" | 生存分组 | 选择 | 给了 | D14 |
| "HSC-like subtype exhibited significantly worse prognosis compared to GMP-like" | 生存 | KM | 给了方向 | — |

---

## ④ 逐步可运行代码

### 4.1 微调

```python
# code/10_01_finetune.py
import torch, torch.nn as nn
# 承接 03/04 的 model 实例（model = ScTransMILCode(...)）
# 1) 载入泛癌预训练权重（作者是 CLAM_SB 的 state_dict，键为 bert_encoder.*/attention_net.*/classifiers.*）
#    若用本项目的 ScTransMILCode（encoder.*/attention.*/cls.*），必须先做 key remap：
#    state = {k.replace("bert_encoder.", "encoder.").replace("attention_net.", "attention.")
#              .replace("classifiers.", "cls."): v for k, v in state.items()}
state = torch.load("raw/pretrained/sctransmil_bicls/model.pt", map_location="cpu")
model.load_state_dict(state, strict=False)
# 2) 冻结策略（D13，与 utils/finetune_utils.py 的 def_freeze_layers=4 一致）：
#    冻结 embeddings + encoder.layer[0:4]；训练 layer[4:6] + pooler + aggregator
for name, p in model.named_parameters():
    p.requires_grad = False
for name, p in model.named_parameters():
    if any(k in name for k in ["encoder.encoder.layer.4", "encoder.encoder.layer.5",
                               "encoder.pooler", "attention", "cls"]):
        p.requires_grad = True
# 3) 二分类微调，30 epochs，早停
```

> 层名以实际 `BertModel` 的 `encoder.layer.{i}` 为准；先 `print(model)` 确认。若权重是作者 `CLAM_SB` 的 state_dict，key 名会不同。

### 4.2 三种策略对比

```python
# code/10_02_strategies.py
strategies = {
    "no_ft": [],
    "agg_only": ["attention", "cls"],
    "last2_agg": ["encoder.encoder.layer.4", "encoder.encoder.layer.5", "attention", "cls"],
}
```

### 4.3 亚型

```python
# code/10_03_subtypes.py
import numpy as np, pandas as pd
from sklearn.neighbors import KNeighborsClassifier
# 输入（由 05_04/微调模型导出）：
#   emb:      (N_cells, 256) 实例嵌入
#   scores:   (N_cells,)     **softmax 后** cell score（D22）
#   cell_type:(N_cells,)     原始细胞类型（含正常分化状态）
#   malignant:(N_cells,)     布尔，来自 AML 原始研究的 malignant 注释（仅用于评估，不参与训练）
NORMAL_STATES = ["HSC", "GMP", "Prog", "ProMono", "Mono", "cDC"]
is_normal = np.isin(cell_type, NORMAL_STATES)
knn = KNeighborsClassifier(n_neighbors=15).fit(emb[is_normal], cell_type[is_normal])
mal = scores > 0.2
subtype = knn.predict(emb[mal])
out = pd.DataFrame({"cell_id": np.arange(len(scores))[mal],
                    "tumor_score": scores[mal], "subtype": subtype})
out.to_csv("results/aml_subtypes.csv", index=False)
print(out.subtype.value_counts())
```

### 4.4 纵向分析

```python
# code/10_04_longitudinal.py
# 按 sample 的 treatment time 排序，统计每个时间点的肿瘤细胞数与亚型比例
```

### 4.5 生存

```python
# code/10_05_survival.py
import pandas as pd, numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

TOP_GENES, TOP_PATIENTS = 30, 50
bulk = pd.read_csv("raw/bulk_aml_dbgap/expr.tsv", sep="\t", index_col=0)  # 行=基因，列=患者
surv = pd.read_csv("raw/bulk_aml_dbgap/survival.tsv", sep="\t")           # patient_id,time,event
att = pd.read_csv("results/attention_aml_subtypes.tsv", sep="\t")        # subtype,gene,attention_mean

picked, used = {}, set()
for st, grp in att.groupby("subtype"):
    genes = grp.sort_values("attention_mean", ascending=False).head(TOP_GENES).gene.tolist()
    genes = [g for g in genes if g in bulk.index]
    score = bulk.loc[genes].mean(axis=0)
    cand = score.sort_values(ascending=False).index.difference(list(used))[:TOP_PATIENTS]
    picked[st] = cand; used.update(cand)

rows = []
for st, pts in picked.items():
    for p in pts:
        r = surv[surv.patient_id == p].iloc[0]
        rows.append(dict(patient_id=p, subtype=st, time=r.time, event=r.event))
df = pd.DataFrame(rows)
df.to_csv("results/survival_aml.csv", index=False)
kmf = KaplanMeierFitter()
for st, sub in df.groupby("subtype"):
    kmf.fit(sub.time, sub.event, label=st); kmf.plot_survival_function()
# HSC-like vs GMP-like
a = df[df.subtype.str.startswith("HSC")]; b = df[df.subtype.str.startswith("GMP")]
print(logrank_test(a.time, b.time, a.event, b.event))
```

### 4.6 敏感性分析

```python
# code/10_06_sensitivity.py
# 变化：top 30/50 基因；top 30/50/100 患者；score 阈值 0.1/0.2/0.3
# 报告 HSC-like vs GMP-like 的 log-rank p 是否稳定
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| 微调策略 | 最后 2 块 + aggregator | 只 aggregator/全量 | 性能与过拟合 |
| 微调 epoch | 30 | 10–50 | 早停 |
| 学习率 | 2e-4 | 1e-4–5e-4 | 收敛 |
| tumor score 阈值 | 0.2 | 0.1–0.5 | 恶性细胞数/亚型 |
| kNN k | 15 | 5–30 | 亚型边界 |
| 生存基因数 | top 30 | 20–100 | 分组稳定性 |
| 每组患者数 | top 50 | 30–100 | 统计功效 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| 数据规模 | 35+6 / 30k / 21 donors | ±5% | ❌ |
| 微调划分 | 20+3 / 15+3 | 必须 | ❌ |
| 最佳策略 | 最后 2 块 + aggregator | 方向一致 | ⚠️ |
| 10+3 平台 | 是 | — | ⚠️ |
| 亚型数 | 6 | 5–7 | ⚠️ |
| 每亚型 >1000 细胞 | 是 | — | ⚠️ |
| HSC-like vs GMP-like | HSC 更差，p<0.05 | 方向一致 | ❌ |
| 生存分组 | 6×50=300 | — | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 微调后全判 tumor | donor setting 标签噪声 | 早停/学习率；报告噪声 |
| 层名不匹配 | state_dict key 不同 | `print(model)` 后按前缀解冻 |
| 亚型数不是 6 | 阈值/kNN k | 敏感性分析；不要硬调 |
| 生存无差异 | 分组不稳定 | 换基因数/患者数；报告 p 的敏感性 |
| dbGaP 无权限 | 未申请 | 标记未复现（D12） |
| 轨迹方向反了 | 正常状态标签顺序 | 按 HSC→Mono 生物学顺序解释 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 5a | `figures/Fig5A_finetune.*` | 三种策略 |
| Fig. 5b | `figures/Fig5BC_aml_embedding.*` | 分化轨迹 |
| Fig. 5c | 同上 | 微调前后 score |
| Fig. 5d | `figures/Fig5DE_subtypes.*` | 亚型比例 |
| Fig. 5e | 同上 | 恶性细胞轨迹 |
| Fig. 5f | `figures/Fig5F_longitudinal.*` | 治疗时间 |
| Fig. S8–S12 | `results/10_*` | 数据/微调/基因/生存 |

> **交付物**：`results/aml_subtypes.csv`、`results/survival_aml.csv`、`figures/Fig5*`。完成后进入 `11_潜在空间_scIB与通路富集.md`。
