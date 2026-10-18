# 11 潜在空间、scIB 与通路富集

> **目标文献**：Tang Z, et al. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* Nature Communications (2026) 17:10800. DOI: 10.1038/s41467-026-77538-5
> **模块对应原文**：Fig. 3g–j；Fig. 4h；Results「… latent space representations …」「… pathway enrichment analysis」；Methods「Visualization」「Functional enrichment analysis」
> **前置依赖**：`05`（bag embedding）、`08`（BRCA 推理）、`09`（attention 基因）
> **复现难度 / 预计工时 / 计算资源**：★★★☆☆ / 2–3 人日 / CPU 为主

---

## ① 科学问题与逻辑位置

- **上游**：训练好的 bag embedding（256 维）、attention top-100 基因。
- **下游**：`figures/Fig3G–J`、`Fig4H`；13 号批判性评述。
- **本模块回答两个问题**：
  1. **表征**：scTransMIL 学到的 bag embedding 是否既保留癌种身份、又混合数据集批次？（scIB 量化）
  2. **通路**：attention top-100 基因富集到哪些通路？是否与疾病生物学一致？

---

## ② 背景原理

### 2.1 scIB 的八个子指标

scIB（Luecken et al. 2022）把整合质量分为两类：

**Bio Conservation（生物保守）**

| 指标 | 含义 | 方向 |
|---|---|---|
| Isolated labels | 稀有群体是否被分开 | 越高越好 |
| KMeans NMI | 无监督聚类与真实标签的互信息 | 越高越好 |
| KMeans ARI | 聚类与标签的调整 Rand 指数 | 越高越好 |
| Silhouette label | 同类样本的紧致/分离 | 越高越好 |

**Batch Correction（批次校正）**

| 指标 | 含义 | 方向 |
|---|---|---|
| Silhouette batch | 批次间混合程度 | 越高越好（接近 0/负） |
| KBET | kNN 邻域内批次分布是否随机 | 越高越好 |
| Graph connectivity | 同类型跨批次是否连通 | 越高越好 |
| PCR comparison | 批次可解释方差的比例 | 越低越好（scIB 取反） |

**Total**：Bio 与 Batch 的加权平均（scIB 默认）。

### 2.2 本项目的"非常规"设定

原文：

> "We use sample embedding as cell embedding, cancer type as cell type, dataset ID (2 datasets: Pan-cancer and BRCA) as batch ID."

也就是说：

- "cell" = sample（bag embedding）；
- "cell type" = cancer type（29/30 类）；
- "batch" = 数据集（仅 2 个：泛癌 vs BRCA）。

**注意**：只有 2 个 batch，KBET/PCR 的统计意义有限；scIB 的分数应理解为"相对排序"，不是绝对质量。

### 2.3 为什么 scTransMIL 的 Bio Conservation 反而较低

原文解释：

> "scTransMIL's instance encoder discards specific gene expression signals … mitigates batch effects … attenuates the biological heterogeneity among samples, resulting in diminished biological conservation."

即：rank 编码丢弃绝对值 → 批次效应↓、但细粒度生物差异也↓。这是一个**明确的 trade-off**，不是缺陷。

### 2.4 通路富集

- 工具：GSEApy `enrichr`；
- 检验：Fisher 单侧；
- 校正：BH；
- 输入：attention top-100 基因（SYMBOL）；
- 输出：top-20 通路，FDR<0.001。

**背景基因集**：enrichr 默认用其数据库的全基因组背景；本项目应确认与词表一致（D15）。若基因是 Ensembl，必须先转 SYMBOL。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "we analyzed the latent space representations … by extracting and visualizing bag-level embeddings" | bag UMAP | scanpy | — | — |
| "scIB … sample embedding as cell embedding, cancer type as cell type, dataset ID … as batch ID" | scIB 设定 | scib | 给了 | — |
| "scTransMIL achieved a higher Batch Correction score and a Total Aggregate score" | Batch/Total 更高 | scIB | 未给数值 | 以 Source Data 为准 |
| "scTransMIL achieves superior batch correction but exhibits poorer biological conservation" | trade-off | 解释 | — | — |
| "Functional enrichment … GSEApy (gseapy.enrichr) … one-sided Fisher's exact test … BH" | 富集 | gseapy | 给了 | D15 |
| "top 20 enriched pathways (FDR < 0.001)" | top-20 | 过滤 | 给了 | — |

---

## ④ 逐步可运行代码

### 4.1 scIB 评估

```python
# code/11_01_scib.py
import os, numpy as np, pandas as pd, scib
from sklearn.cluster import KMeans

# adata: n_samples × 256（bag embedding），obs 含 cancer_type / dataset
adata.obs["kmeans"] = KMeans(n_clusters=adata.obs.cancer_type.nunique(),
                             random_state=20261018, n_init=10).fit_predict(adata.X)

res = {}
# Bio conservation（scIB 的函数名以安装版本为准）
res["isolated_labels"] = scib.metrics.isolated_labels(
    adata, label_key="cancer_type", batch_key="dataset", embed="X")
res["nmi"] = scib.metrics.nmi(adata, adata, cluster_key="kmeans", label_key="cancer_type")
res["ari"] = scib.metrics.ari(adata, adata, cluster_key="kmeans", label_key="cancer_type")
res["silhouette_label"] = scib.metrics.silhouette(adata, label_key="cancer_type", embed="X")
# Batch correction
res["silhouette_batch"] = scib.metrics.silhouette_batch(
    adata, batch_key="dataset", label_key="cancer_type", embed="X")
for name, fn in [("kbet", scib.metrics.kBET),
                 ("graph_connectivity", scib.metrics.graph_connectivity)]:
    try:
        res[name] = fn(adata, batch_key="dataset", label_key="cancer_type", embed="X") \
            if name == "kbet" else fn(adata, label_key="cancer_type")
    except Exception as e:
        res[name] = np.nan; print("skip", name, e)
pd.Series(res).to_csv("results/11_scib_submetrics.csv")
```

> scIB 的 API 随版本变化（`scib.metrics.metrics` vs `scib.me.metrics`）；以安装版本为准。只有 2 个 batch 时 KBET 可能报错/不稳定，记录并跳过。

### 4.2 UMAP

```python
# code/11_02_umap.py
import scanpy as sc
ad = sc.AnnData(emb)              # n_samples × 256
ad.obs = meta
sc.pp.neighbors(ad, n_neighbors=15, use_rep="X")
sc.tl.umap(ad, random_state=20261018)
sc.pl.umap(ad, color=["cancer_type", "dataset"], save="bag_umap.png")
```

### 4.3 富集

```python
# code/11_03_enrichment.py
import pandas as pd, gseapy as gp
att = pd.read_csv("results/attention_all_tumorB.tsv", sep="\t")
genes = att.sort_values("attention_mean", ascending=False).head(100).gene.dropna().unique().tolist()
enr = gp.enrichr(gene_list=genes,
                 gene_sets=["GO_Biological_Process_2023", "MSigDB_Hallmark_2020"],
                 organism="human", outdir=None)
res = enr.results
res["FDR"] = res["Adjusted P-value"]
top = res[res.FDR < 0.001].sort_values("FDR").head(20)
top.to_csv("results/11_enrichment_top20.csv", index=False)
print(top[["Term", "FDR", "Overlap"]])
```

### 4.4 与 DEA 基因的交集

```python
# code/11_04_overlap.py
att_genes = set(att.head(100).gene)
dea_path = "results/09_dea_tumorB_vs_normalB.csv"
if not os.path.exists(dea_path):
    raise FileNotFoundError("先运行 09_05_dea_compare.py 生成 DEA 结果（契约：results/09_dea_tumorB_vs_normalB.csv）")
dea_genes = set(pd.read_csv(dea_path).gene)
print("overlap:", len(att_genes & dea_genes))
print("BLK in attention:", "BLK" in att_genes, " in DEA:", "BLK" in dea_genes)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| `n_neighbors`（UMAP） | 15 | 5–50 | 视觉 |
| scIB label/batch | cancer_type / dataset | — | 指标定义 |
| 富集库 | GO BP + Hallmark | KEGG/Reactome | 通路 |
| FDR 阈值 | 0.001 | 0.05 | 通路数 |
| 背景 | enrichr 默认 | 词表基因 | p 值 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| bag UMAP 按癌种分簇 | 是 | 定性 | ⚠️ |
| 两个数据集混合 | 部分混合 | 定性 | ⚠️ |
| Batch Correction | scTransMIL 最高 | 相对 | ⚠️ |
| Bio Conservation | 低于 PCA/FM 变体 | 相对 | ⚠️ |
| Total | 最高 | 相对 | ⚠️ |
| top-20 通路 FDR<0.001 | 是 | 方向一致 | ⚠️ |
| 通路与 ALL 相关 | 免疫/细胞因子/抗原受体 | 定性 | ⚠️ |
| BLK 在 attention、不在 DEA | 是 | 定性 | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| scIB 报 "batch too small" | 只有 2 个 batch | 跳过 KBET/PCR，记录 |
| UMAP 按数据集分 | 批次效应 | 正是要量化的；与原文对照 |
| 富集无结果 | 基因是 Ensembl/大小写 | 转 SYMBOL 大写 |
| 背景不一致 | enrichr 默认全基因组 | 记录；必要时自定义背景 |
| Total 分数算错 | Bio/Batch 权重 | 用 scIB 默认 |
| 与原文排序不符 | 数据版本/训练 | 以相对排序为主 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 3g | `figures/Fig3GH_bag_umap.*` | 癌种 |
| Fig. 3h | 同上 | 数据集 |
| Fig. 3i | `results/11_scib_submetrics.csv` | 八项 |
| Fig. 3j | `results/11_scib_total.csv` | Batch/Bio/Total |
| Fig. 4h | `figures/Fig4H_enrichment.*` | top-20 |
| Methods「Visualization/Enrichment」 | 本文件 §2 | — |

> **交付物**：`results/11_scib_*.csv`、`results/11_enrichment_top20.csv`、`figures/Fig3G–J`、`Fig4H`。完成后进入 `12_计算性能与可复现性.md`。
