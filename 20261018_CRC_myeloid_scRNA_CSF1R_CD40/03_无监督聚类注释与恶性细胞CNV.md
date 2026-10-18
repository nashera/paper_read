# 03 无监督聚类、细胞注释与恶性细胞 CNV

> **目标文献**：Zhang L, Li Z, Skrzypczynska KM, et al. *Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer.* Cell. 2020;181(2):442–459.e29. DOI: [10.1016/j.cell.2020.03.048](https://doi.org/10.1016/j.cell.2020.03.048)
> **模块对应原文**：Results p.444–445（Fig. 1A–1C）；Fig. S2、S3；STAR Methods「Unsupervised clustering analysis」（p.e9）、「Dimensionality reduction using t-SNE and UMAP」（p.e9）、「Identification of signature genes」（p.e9）、「Tissue distribution of clusters」（p.e9）、「Similarity analysis of clusters from different datasets」（p.e9–e10）、「Identification of malignant cells with CNV estimation」（p.e10）
> **前置依赖**：`01_环境搭建与数据获取.md`；`02_人类CRC数据接入与统一预处理.md`（`processed/human_10x_norm.rds`、`human_smartseq2_norm.rds`）；`PROJECT_SPEC.md` §4 契约（`human_annot.rds`、`human_cnv.rds`）、§5 术语、§8 D7–D9
> **复现难度 / 预计工时 / 计算资源**：★★★★☆ / 3–5 人日 / 64 GB 内存、16 核 CPU；CNV 可按患者并行

---

## ① 科学问题与逻辑位置

### 1.1 本模块要解决什么

把 02 号输出的归一化矩阵，变成**带细胞类型标签、带组织分布、带恶性判定**的注释对象：

- 每个数据集独立聚类（10x 与 Smart-seq2 不合并）；
- 第一轮 8 大类（T、ILC、B、髓系、上皮、内皮、成纤维、恶性）；
- 第二轮在各大大类内再聚类，得到 38/36 个免疫簇、13 个髓系簇（并集；Smart-seq2 只分辨出 9 个）、4 个 ILC、18 个 T、5 个 B、12 个非免疫簇；
  - 口径：10x 38 = 13 髓系 + 3 ILC + 17 T + 5 B；Smart-seq2 36 = 9 髓系 + 4 ILC + 18 T + 5 B。不要把 `13+4+18+5=40` 当作单平台簇数。
- 用 marker 表把簇映射到原文标签（`hM01–hM13`、`hT01–hT18`、`hB01–hB05`、`hI01–hI04`、`hE01–hE06`、`hF01–hF02`、`hC01–hC04`）；
- 用跨平台 logistic regression 做相似性映射（Fig. S2E）；
- 用 `Ro/e` 计算组织富集（Fig. 2C、S3C）；
- 用 CNV 推断恶性细胞（Fig. S3E–S3F）。

上游：`02` 号的 `human_10x_norm.rds`、`human_smartseq2_norm.rds`。
下游：`04_髓系亚群与TAM二分法.md` 读取 `human_annot.rds` 的髓系子集；`05`–`07` 号模块读取注释与 CNV；14 号用 Q09–Q15 核对。

### 1.2 一条不可违背的纪律

**注释标签一旦确定，下游模块不得重新聚类、不得改标签。** 下游只允许：

- 取子集（如 `cluster_id %in% hM01:hM13`）；
- 在同一标签体系内做差异表达/通路/轨迹；
- 若发现标签错误，回到本模块改，改完更新 `PROJECT_SPEC.md` §4 契约与 `logs/DECISIONS.md`。

原因：04–07 号的所有结果都依赖同一套 `cluster_id`；每个模块各自重聚类会让细胞类型、比例、网络节点全部漂移。

---

## ② 背景原理

### 2.1 图聚类（Seurat）

**步骤**：

1. **HVG**：在 mean expression vs dispersion（variance/mean）上取高变基因；原文约 2,000 个。
2. **PCA**：对 HVG 标准化后做 PCA；取前 10–25 个 PC。
3. **kNN/SNN 图**：每个细胞连到最近的 k 个细胞；SNN 用共享邻居加权。
4. **Louvain/Leiden**：最大化模块度
   $$
   Q=\frac{1}{2m}\sum_{ij}\left[A_{ij}-\frac{k_i k_j}{2m}\right]\delta(c_i,c_j)
   $$
   其中 `A` 是邻接矩阵，`k_i` 是节点度，`m` 是边数，`c_i` 是社区。`resolution` 越大，社区越多（簇越细）。
5. **第一轮/第二轮**：第一轮分辨大类；第二轮在子集内用同样参数再聚类。

**原文参数**：Seurat v2.3.0；约 2,000 HVG；10–25 PCs；resolution 0.6–2.0。复现用 Seurat 5.1.0，保持同样的参数范围，但记录 API 差异（D2/D7）。

### 2.2 跨平台相似性（logistic regression）

原文（Young et al. 2018）：

1. 训练集每个 cluster 下采样到最小簇大小；若最小簇 <50，则下采样到 50。
2. 对每个 cluster `k`，训练二分类 logistic regression：
   $$
   \log\frac{p(y=k\mid x)}{1-p(y=k\mid x)}=\beta_0+\beta^T x
   $$
   用 elastic net（`alpha=0.99`，10-fold CV）。
3. offset：`log(f/(1−f))`，`f` 是该簇在训练集中的比例；测试时 offset=0。
4. 对测试集每个细胞，每个模型输出 logit；在测试 cluster 内平均 → 相似性矩阵；转概率可视化。

**直觉**：每个 cluster 一个"分类器"；测试 cluster 的细胞平均更像哪个训练 cluster，就映射到哪个。

### 2.3 t-SNE / UMAP

- **t-SNE**：在高维空间用高斯核定义相似度 `p_{j|i}`，在低维空间用 t 分布定义 `q_{ij}`，最小化 `KL(P||Q)`。强调局部结构，全局距离不可解释。
- **UMAP**：基于黎曼几何与 fuzzy topological representation，最小化交叉熵；`n_neighbors` 控制局部/全局平衡，`min_dist` 控制点间距。
- 原文：`RunTSNE`/`RunUMAP`，PCs 与聚类相同；UMAP `spread` 1–2。
- **CCA 整合只用于可视化**（Fig. 1C）：`IntegrateData` 后的矩阵不用于聚类或下游分析。

### 2.4 Ro/e（组织富集）

对 cluster × tissue 的列联表，卡方检验给出期望数：
$$
E_{ct}=\frac{O_{c\cdot}O_{\cdot t}}{O_{\cdot\cdot}},\qquad
\mathrm{Ro/e}_{ct}=\frac{O_{ct}}{E_{ct}}
$$
`Ro/e>1` 表示该 cluster 在该组织富集。原文用 `Ro/e` 判定 tumor/normal/blood 富集；Fig. S8B 用 `+++`（>2）、`++`（1–2）、`+`（0.5–1）、`+/−`（0–0.5）、`−`（0）。

### 2.5 CNV 推断（恶性细胞）

原文公式（完整复现）：

1. 过滤低表达基因：`log2(mean(TPM_i,·)+1) < 3.5` 丢弃。
2. `E_ij = log2(TPM_ij/10 + 1)`。
3. 基因中心化：`Er_ij = E_ij − mean_i(E_i,·)`。
4. 截断：`Er ∈ [−3,3]`。
5. 按染色体位置排序，滑动窗口 100 个基因取平均 → `CNV_i`。
6. 每个细胞中心化：`CNV_c = CNV_i − mean_j(CNV_i,j)`（每个细胞均值为 0）。
7. 用非恶性上皮簇 `hE02–hE06` 的 `CNV_c` 均值作为 baseline；每个窗口取 `BaseMax`、`BaseMin`：
   $$
   \mathrm{CNV}_f(i,j)=
   \begin{cases}
   \mathrm{CNV}_c(i,j)-\mathrm{BaseMax}(j), & \mathrm{CNV}_c(i,j)>\mathrm{BaseMax}(j)+0.2\\
   \mathrm{CNV}_c(i,j)-\mathrm{BaseMin}(j), & \mathrm{CNV}_c(i,j)<\mathrm{BaseMin}(j)-0.2\\
   0, & \text{otherwise}
   \end{cases}
   $$
8. 画图时 `CNV_f ∈ [−1,1]`，`[−0.1,0.1]` 置 0。
9. `CNV score = mean(CNV_f²)`；`CNV correlation = Pearson(CNV_f(cell), mean CNV_f(malignant cells))`。

**直觉**：恶性细胞有大片段扩增/缺失；用非恶性上皮做 baseline 去掉组织/平台背景。`CNV score` 高 + `CNV correlation` 高 → 恶性。`infercnv` 是另一套实现，可作敏感性分析（D9）。

---

## ③ 原文逐句拆解

| 原文句（节选） | 中文转述 | 我的操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "we used an unsupervised graph-based clustering algorithm implemented in Seurat v2 (version 2.3.0) for each dataset separately." | 每个数据集独立聚类 | 10x/SS2 分别跑 | 给了版本 | D2/D7 |
| "highly variable genes ... appropriate threshold of mean expression and dispersion (variance/mean)." | HVG 用 mean/dispersion | `FindVariableFeatures(selection.method="mean.var.plot")` | 文献未报告 | 约 2,000 HVG（原文下一句） |
| "PCA was performed on about 2000 variable genes." | 2,000 HVG | `nfeatures=2000` | 给了数量 | 照做 |
| "FindClusters on 10-25 PCs with resolution 0.6-2." | PCs/resolution 范围 | 扫描 resolution，锚定簇数 | 给了范围 | D7：以 13/15 簇为锚 |
| "eight major cell types ... 4 immune ... 4 non-immune ..." | 8 大类 | 第一轮聚类 + marker | 给了类别 | 照做 |
| "malignant cell clusters were also characterized by the exceptionally high amounts of expressed genes." | 恶性细胞基因数高 | 辅助判断 | 文献未报告 | 用 CNV + 上皮 marker |
| "38 and 36 leukocyte cell clusters" | 38/36 簇 | 第二轮聚类 | 给了簇数 | Q09/Q10 |
| "13 myeloid cell clusters, 4 ILC clusters, 18 T cell clusters, and 5 B cell clusters" | 13/4/18/5 | 注释 | 给了簇数 | Q11–Q12 |
| "12 cell clusters"（non-immune） | 12 非免疫簇 | 注释 | 给了簇数 | Q13 |
| "Ro/e = Observed/Expected ... expected ... Chi-square test ... enriched if Ro/e>1" | Ro/e | `03_04_ro_e.R` | 给了公式/阈值 | 照做 |
| "logistic regression ... elastic net ... alpha=0.99 ... 10-fold CV ... offset log(f/(1−f))" | 相似性 | `03_03_similarity.R` | 给了参数 | 照做 |
| "downsampled each cluster to the minimum size ... clusters lower than 50 ... set as 50" | 下采样 | `downsample_clusters()` | 给了规则 | 照做 |
| "CNV ... moving average ... sliding window of 100 genes ... baseline hE02–hE06 ... ±0.2" | CNV | `03_05_cnv.R` | 给了公式/参数 | 照做 |
| "CNV score ... mean squares ... CNV correlation ... Pearson" | 恶性判定 | `03_05_cnv.R` | 给了定义 | 照做 |

---

## ④ 逐步可运行代码

### 4.1 聚类（10x 与 Smart-seq2 分别）

```r
#!/usr/bin/env Rscript
# code/03_01_cluster.R
# 输入: processed/human_10x_norm.rds, processed/human_smartseq2_norm.rds
# 输出: processed/human_10x_clustered.rds, processed/human_smartseq2_clustered.rds

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(SummarizedExperiment)
  library(Matrix); library(readr); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

cluster_one <- function(sce_path, out_path, tag, dims = 1:20) {
  sce <- readRDS(sce_path)
  logtpm <- as.matrix(assay(sce, "logtpm"))
  counts <- as.matrix(assay(sce, "counts"))

  # Seurat 5 对象：counts 用于对象；data 用 logTPM
  obj <- CreateSeuratObject(counts = counts, project = tag,
                            min.cells = 0, min.features = 0)
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = logtpm)

  # HVG：原文 mean/dispersion；Seurat 5 用 "mean.var.plot"
  obj <- FindVariableFeatures(obj, selection.method = "mean.var.plot",
                              nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30,
                verbose = FALSE, seed.use = 20261018)

  # 第一轮：大类（resolution 0.6）
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE,
                      random.seed = 20261018)
  obj$cluster_first <- Idents(obj)

  # 扫描 resolution，记录簇数（供 D7 选择）
  res_seq <- seq(0.6, 2.0, by = 0.2)
  cluster_counts <- sapply(res_seq, function(r) {
    o <- FindClusters(obj, resolution = r, verbose = FALSE,
                      random.seed = 20261018)
    length(levels(Idents(o)))
  })
  names(cluster_counts) <- res_seq
  write_tsv(data.frame(resolution = res_seq, n_clusters = cluster_counts),
            file.path(proj, "results", paste0("03_cluster_counts_", tag, ".tsv")))

  # 主分析取 resolution=1.0（可按簇数微调；记录 D7）
  obj <- FindClusters(obj, resolution = 1.0, verbose = FALSE,
                      random.seed = 20261018)
  obj$cluster_id <- Idents(obj)

  # 嵌入
  obj <- RunUMAP(obj, dims = dims, seed.use = 20261018, spread = 1)
  obj <- RunTSNE(obj, dims = dims, seed.use = 20261018)

  saveRDS(obj, out_path)
  cat("[OK]", tag, ":", ncol(obj), "cells; first-round clusters:",
      length(unique(obj$cluster_first)), "; resolution 1.0 clusters:",
      length(unique(obj$cluster_id)), "\n")
  obj
}

obj10 <- cluster_one(file.path(proj, "processed", "human_10x_norm.rds"),
                     file.path(proj, "processed", "human_10x_clustered.rds"), "10x")
objss <- cluster_one(file.path(proj, "processed", "human_smartseq2_norm.rds"),
                     file.path(proj, "processed", "human_smartseq2_clustered.rds"), "ss2")
```

> ⚠️ **导师批注**：Seurat 5 的 `FindClusters` 默认算法是 Louvain（`algorithm=1`），与 Seurat v2 的实现不完全相同；若要更接近原文，可设 `algorithm=1, method="igraph"`，但不同 igraph 版本仍可能有差异。**不要**为了凑 38/36 簇反复改 resolution 直到"看起来对"；应先用 marker 稳定性选分辨率，再接受实际簇数（L0 §12）。

### 4.2 注释：marker 表与标签映射

```r
#!/usr/bin/env Rscript
# code/03_02_annotate.R
# 输入: processed/human_*_clustered.rds
# 输出: processed/human_annot.rds（10x 主注释；Smart-seq2 标签映射见 §4.3）
#        results/03_cluster_markers.tsv, results/03_cluster_labels.tsv

suppressPackageStartupMessages({
  library(Seurat); library(readr); library(dplyr); library(tidyr)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

# ---- marker 表（原文 Fig. 1–2、S2–S4） ----
markers <- tribble(
  ~major, ~label, ~genes,
  "Myeloid", "hM01_Mast-TPSAB1",       "TPSAB1,TPSB2,CPA3,MS4A2,KIT,HDC",
  "Myeloid", "hM02_pDC-LILRA4",        "LILRA4,LILRB4,IL3RA,GZMB,TCF4,IRF7",
  "Myeloid", "hM03_cDC2-CD1C",         "CD1C,FCER1A,CLEC10A,CD1E,CD1A",
  "Myeloid", "hM04_cDC1-BATF3",        "BATF3,XCR1,CLEC9A,IRF8,IDO1",
  "Myeloid", "hM05_Mono-CD14",         "CD14,S100A8,S100A9,VCAN,LYZ",
  "Myeloid", "hM06_Mono-CD16",         "FCGR3A,MS4A7,CX3CR1,LILRB1",
  "Myeloid", "hM07_Mono-CD14CD16",     "CD14,FCGR3A,MS4A7,IFITM2",
  "Myeloid", "hM08_Macro-NLRP3",       "NLRP3,IL1B,EREG,ITGAX",
  "Myeloid", "hM09_Macro-PLTP",        "PLTP,LYVE1,IL10,APOE",
  "Myeloid", "hM10_Macro-IL1B",        "IL1B,CCL3,VEGFA,APOE",
  "Myeloid", "hM11_Monolike-FCN1",     "FCN1,CD14,S100A8,C1QC,SPP1",
  "Myeloid", "hM12_TAM-C1QC",          "C1QA,C1QB,C1QC,TREM2,MERTK,CD80",
  "Myeloid", "hM13_TAM-SPP1",          "SPP1,MARCO,VEGFA,FN1,PPARG",
  "ILC",     "hI01_NK-CD16",           "FCGR3A,GNLY,NKG7,GZMB",
  "ILC",     "hI02_NK-GZMK",           "GZMK,GZMA,NKG7,IL7R",
  "ILC",     "hI03_NK-CD103",          "ITGAE,CD160,LAYN,GNLY",
  "ILC",     "hI04_ILC3-SLC4A10",      "SLC4A10,IL23R,RORC,KIT",
  "T",       "hT01_CD4-CCR7",          "CD4,CCR7,LEF1,SELL",
  "T",       "hT02_CD4-ANXA1",         "CD4,ANXA1,IL7R,ITM2A",
  "T",       "hT03_CD4-GNLY",          "CD4,GNLY,NKG7,GZMA",
  "T",       "hT04_CD4-TCF7",          "CD4,TCF7,LEF1,SELL",
  "T",       "hT05_CD4-CXCR6",         "CD4,CXCR6,ITGAE,CD69",
  "T",       "hT06_CD4-CXCR5",         "CD4,CXCR5,BCL6,PDCD1",
  "T",       "hT07_CD4-GZMK",          "CD4,GZMK,CD27,IL7R",
  "T",       "hT08_CD4-IL23R",         "CD4,IL23R,CCR6,RORC",
  "T",       "hT09_CD4-CXCL13",        "CD4,CXCL13,PDCD1,CTLA4",
  "T",       "hT10_CD4-FOXP3",         "CD4,FOXP3,IL2RA,CTLA4",
  "T",       "hT11_CD4-CTLA4",         "CD4,CTLA4,PDCD1,TIGIT",
  "T",       "hT12_CD8-LEF1",          "CD8A,LEF1,CCR7,SELL",
  "T",       "hT13_CD8-GPR183",        "CD8A,GPR183,ITM2A,CD69",
  "T",       "hT14_CD8-CX3CR1",        "CD8A,CX3CR1,FCGR3A,GZMB",
  "T",       "hT15_CD8-GZMK",          "CD8A,GZMK,GZMA,CD27",
  "T",       "hT16_CD8-CD6",           "CD8A,CD6,IL7R,CCR7",
  "T",       "hT17_CD8-CD160",         "CD8A,CD160,LAYN,TIGIT",
  "T",       "hT18_CD8-LAYN",          "CD8A,LAYN,PDCD1,CTLA4",
  "B",       "hB01_PlasmaB-IgG",       "MZB1,JCHAIN,IGHG1,IGHG3",
  "B",       "hB02_GALTB-IgA",         "GALTB,IGHA1,IGHA2,JCHAIN",
  "B",       "hB03_FollicularB-IgD",   "IGHD,MS4A1,CD79A,FCER2",
  "B",       "hB04_FollicularB-MS4A1", "MS4A1,CD79A,CD19,CD37",
  "B",       "hB05_GCBCell-LRMP",      "LRMP,TCL1A,BANK1,FCER2",
  "Endothelial", "hE01_Endothelium-ACKR1", "ACKR1,CDH5,ERG,PLVAP",
  "Epithelial", "hE02_Enterocyte-GUCA2B", "GUCA2B,SLC26A3,CDX2,KRT20",
  "Epithelial", "hE03_Goblet-SPINK4",  "SPINK4,TFF3,SPDEF,MUC2",
  "Epithelial", "hE04_Goblet-REG4",    "REG4,CLCA4,MUC2,TFF3",
  "Epithelial", "hE05_Stemlike-LGR5",  "LGR5,ASCL2,OLFM4,SMOC2",
  "Epithelial", "hE06_UnIdent",        "EPCAM,KRT8,KRT18",
  "Fibroblast", "hF01_Myofib-ACTA2",   "ACTA2,TAGLN,MYH11,COL1A1",
  "Fibroblast", "hF02_CAF-FAP",        "FAP,DCN,COL1A1,COL3A1,MMP2"
)
write_tsv(markers, file.path(proj, "results", "03_marker_panel.tsv"))

# ---- 计算每个 cluster 的 marker 平均表达与模块分数 ----
annotate_obj <- function(obj_path, out_path, tag) {
  obj <- readRDS(obj_path)
  logtpm <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
  clusters <- obj$cluster_id

  score_mat <- matrix(NA_real_, nrow = length(unique(clusters)),
                      ncol = nrow(markers),
                      dimnames = list(sort(unique(clusters)), markers$label))
  for (i in seq_len(nrow(markers))) {
    g <- intersect(strsplit(markers$genes[i], ",")[[1]], rownames(logtpm))
    if (length(g) == 0) next
    score_mat[, i] <- tapply(colMeans(logtpm[g, , drop = FALSE]), clusters, mean)
  }

  # 每个 cluster 取最高分标签（仅作初稿；最终需人工核对）
  best <- apply(score_mat, 1, function(x) {
    if (all(is.na(x))) return(NA_character_)
    colnames(score_mat)[which.max(x)]
  })
  label_df <- data.frame(cluster_id = rownames(score_mat), best_label = best)
  write_tsv(label_df, file.path(proj, "results", paste0("03_cluster_labels_", tag, ".tsv")))

  # 写入对象（初稿标签；人工复核后覆盖）
  obj$major_type <- markers$major[match(best, markers$label)]
  obj$cluster_label_draft <- best
  saveRDS(obj, out_path)
  cat("[OK]", tag, "注释初稿完成；请人工核对 results/03_cluster_labels_", tag, ".tsv\n", sep = "")
}

annotate_obj(file.path(proj, "processed", "human_10x_clustered.rds"),
             file.path(proj, "processed", "human_10x_annot_draft.rds"), "10x")
annotate_obj(file.path(proj, "processed", "human_smartseq2_clustered.rds"),
             file.path(proj, "processed", "human_smartseq2_annot_draft.rds"), "ss2")
```

> ⚠️ **导师批注**：自动取最高分只是**初稿**。原文的标签是人工根据 marker、组织分布、跨平台相似性共同确定的；必须逐簇核对 `03_cluster_markers.tsv`、`03_cluster_labels_*.tsv` 与 Fig. 2/S2/S4 的 marker 热图。最终标签手工写入 `results/03_cluster_labels_final.tsv`（10x）与 `results/03_smartseq2_cluster_labels_final.tsv`（Smart-seq2），再由 `03_02b_finalize_annot.R` 覆盖对象。不要用自动标签直接跑下游。

```r
#!/usr/bin/env Rscript
# code/03_02b_finalize_annot.R
# 输入: processed/human_10x_annot_draft.rds + results/03_cluster_labels_final.tsv
#       processed/human_smartseq2_annot_draft.rds + results/03_smartseq2_cluster_labels_final.tsv
# 输出: processed/human_annot.rds（10x，SCE）, processed/human_smartseq2_annot.rds（Smart-seq2，SCE）
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")

finalize_one <- function(draft_path, final_path, out_path, platform) {
  obj <- readRDS(draft_path)
  final <- read_tsv(final_path, show_col_types = FALSE)
  stopifnot(all(c("cluster_id", "final_label", "major_type") %in% colnames(final)))
  m  <- setNames(final$final_label, as.character(final$cluster_id))
  mt <- setNames(final$major_type,  as.character(final$cluster_id))
  # draft 对象里保存的是 Seurat 聚类编号；不要用尚未存在的 obj$cluster_id
  draft_ids <- if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    as.character(obj$seurat_clusters)
  } else {
    as.character(obj$cluster_label_draft)
  }
  obj$cluster_id <- unname(m[draft_ids])
  obj$major_type <- unname(mt[draft_ids])
  obj$platform  <- platform
  if (!"prior_T_cell" %in% colnames(obj@meta.data)) obj$prior_T_cell <- FALSE
  stopifnot(!any(is.na(obj$cluster_id)), !any(is.na(obj$major_type)))
  sce <- as.SingleCellExperiment(obj)
  stopifnot(inherits(sce, "SingleCellExperiment"))
  # 统一提供 logtpm assay：Seurat 的 data 已是 log-normalized；若没有 data 则从 counts 计算
  if (!"logtpm" %in% assayNames(sce)) {
    if ("data" %in% assayNames(sce)) {
      assay(sce, "logtpm") <- assay(sce, "data")
    } else if ("logcounts" %in% assayNames(sce)) {
      assay(sce, "logtpm") <- assay(sce, "logcounts")
    } else {
      assay(sce, "logtpm") <- log2(t(t(assay(sce, "counts")) /
        colSums(assay(sce, "counts"))) * 1e6 + 1)
    }
  }
  if (!"tpm" %in% assayNames(sce) && "logtpm" %in% assayNames(sce)) {
    assay(sce, "tpm") <- pmax(2^assay(sce, "logtpm") - 1, 0)
  }
  saveRDS(sce, out_path)
  cat("[OK]", platform, "最终注释写入", basename(out_path), ":", ncol(sce), "cells\n")
}

finalize_one(
  file.path(proj, "processed", "human_10x_annot_draft.rds"),
  file.path(proj, "results",   "03_cluster_labels_final.tsv"),
  file.path(proj, "processed", "human_annot.rds"),
  "10x"
)
finalize_one(
  file.path(proj, "processed", "human_smartseq2_annot_draft.rds"),
  file.path(proj, "results",   "03_smartseq2_cluster_labels_final.tsv"),
  file.path(proj, "processed", "human_smartseq2_annot.rds"),
  "Smart-seq2"
)
```

### 4.3 跨平台相似性（logistic regression）

```r
#!/usr/bin/env Rscript
# code/03_03_similarity.R
# 输入: processed/human_10x_annot_draft.rds（训练）, human_smartseq2_annot_draft.rds（测试）
# 输出: results/03_cross_platform_similarity.tsv, figures/FigS2E_similarity.pdf/.png

suppressPackageStartupMessages({
  library(Seurat); library(glmnet); library(readr); library(dplyr); library(pheatmap)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

train <- readRDS(file.path(proj, "processed", "human_10x_annot_draft.rds"))
test  <- readRDS(file.path(proj, "processed", "human_smartseq2_annot_draft.rds"))
Xtr <- t(as.matrix(GetAssayData(train, assay = "RNA", layer = "data")))
Xte <- t(as.matrix(GetAssayData(test,  assay = "RNA", layer = "data")))
common <- intersect(colnames(Xtr), colnames(Xte))
Xtr <- Xtr[, common, drop = FALSE]; Xte <- Xte[, common, drop = FALSE]

# 下采样：每簇到最小簇大小（下限 50；不足 50 的簇全用，不做有放回过采样）
downsample_clusters <- function(X, y, min_size = 50) {
  idx <- split(seq_len(nrow(X)), y)
  target <- max(min_size, min(sapply(idx, length)))
  out <- unlist(lapply(idx, function(i) sample(i, min(length(i), target))))
  list(X = X[out, , drop = FALSE], y = y[out])
}

ytr <- factor(train$cluster_id)
ds <- downsample_clusters(Xtr, ytr)
Xtr <- ds$X; ytr <- ds$y

# 每个 cluster 一个二分类器
fits <- list()
for (k in levels(ytr)) {
  y <- as.integer(ytr == k)
  f <- mean(y)
  off <- rep(log(f / (1 - f)), length(y))
  cv <- cv.glmnet(Xtr, y, family = "binomial", alpha = 0.99,
                  offset = off, nfolds = 10, standardize = TRUE)
  fits[[k]] <- cv
}

# 测试集 logit（offset=0），按测试 cluster 平均
logit_mat <- sapply(fits, function(fit) {
  as.numeric(predict(fit, newx = Xte, type = "link", s = "lambda.min",
                     offset = rep(0, nrow(Xte))))
})
colnames(logit_mat) <- names(fits)
test_cluster <- test$cluster_id
sim <- rowsum(logit_mat, group = test_cluster) /
  as.numeric(table(test_cluster))
prob <- 1 / (1 + exp(-sim))

write_tsv(as.data.frame(prob) %>% tibble::rownames_to_column("test_cluster"),
          file.path(proj, "results", "03_cross_platform_similarity.tsv"))

pdf(file.path(proj, "figures", "FigS2E_similarity.pdf"), width = 8, height = 7)
pheatmap(prob, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "10x (columns) vs Smart-seq2 (rows) cluster similarity")
dev.off()
png(file.path(proj, "figures", "FigS2E_similarity.png"), width = 8, height = 7,
    units = "in", res = 300)
pheatmap(prob, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "10x (columns) vs Smart-seq2 (rows) cluster similarity")
dev.off()
```

> ⚠️ **导师批注**：`cv.glmnet` 的 `offset` 参数在 glmnet 4.x 中可用；若版本不支持，改用 `glmnet(..., offset=)` + 手动 CV。`predict` 的 `s="lambda.min"` 与原文的 CV 选择一致。相似性矩阵的**行/列方向**要写清：本代码行=Smart-seq2 测试簇、列=10x 训练簇；不要在图注里写反。

### 4.4 Ro/e

```r
#!/usr/bin/env Rscript
# code/03_04_ro_e.R
# 输入: processed/human_annot.rds（含 cluster_id, tissue）
# 输出: results/ro_e.tsv, figures/Fig2C_ro_e.pdf/.png

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
# 两个平台分别计算 Ro/e，再合并；06 号的“双平台一致”过滤依赖 platform 列
compute_roe <- function(path, platform) {
  obj <- readRDS(path)
  stopifnot(inherits(obj, "SingleCellExperiment"))
  meta <- as.data.frame(colData(obj))
  stopifnot(all(c("cluster_id", "tissue") %in% colnames(meta)))
  tab <- table(meta$cluster_id, meta$tissue)
  chi <- suppressWarnings(chisq.test(tab))
  expected <- chi$expected
  as.data.frame(tab) %>%
    rename(cluster_id = Var1, tissue = Var2, observed = Freq) %>%
    mutate(expected = as.numeric(expected),
           ro_e = observed / expected,
           platform = platform)
}
ro_e <- bind_rows(
  compute_roe(file.path(proj, "processed", "human_annot.rds"), "10x"),
  compute_roe(file.path(proj, "processed", "human_smartseq2_annot.rds"), "Smart-seq2")
)
write_tsv(ro_e, file.path(proj, "results", "ro_e.tsv"))

p <- ggplot(ro_e, aes(tissue, cluster_id, fill = ro_e)) +
  geom_tile() + scale_fill_gradient2(midpoint = 1, low = "blue", high = "red") +
  labs(title = "Ro/e (observed/expected)", x = NULL, y = NULL) +
  theme_minimal(base_size = 8)
ggsave(file.path(proj, "figures", "Fig2C_ro_e.pdf"), p, width = 6, height = 10)
ggsave(file.path(proj, "figures", "Fig2C_ro_e.png"), p, width = 6, height = 10, dpi = 300)
```

### 4.5 CNV

```r
#!/usr/bin/env Rscript
# code/03_05_cnv.R
# 输入: processed/human_smartseq2_annot.rds（Smart-seq2 CD45- 非免疫细胞）
# 输出: processed/human_cnv.rds, figures/FigS3E_cnv.pdf/.png, results/03_cnv_scores.tsv

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(GenomicFeatures)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(Matrix); library(readr); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

sce <- readRDS(file.path(proj, "processed", "human_smartseq2_annot.rds"))
# 只用 CD45- 非免疫细胞（major_type 为 Epithelial/Endothelial/Fibroblast/Malignant）
sce <- sce[, sce$major_type %in% c("Epithelial", "Fibroblast", "Malignant")]
tpm <- as.matrix(assay(sce, "tpm"))

# 1. 过滤低表达基因
keep <- log2(rowMeans(tpm) + 1) >= 3.5
tpm <- tpm[keep, , drop = FALSE]

# 2. E = log2(TPM/10 + 1)
E <- log2(tpm / 10 + 1)

# 3. 基因中心化
Er <- E - rowMeans(E)

# 4. 截断 [-3,3]
Er[Er > 3] <- 3; Er[Er < -3] <- -3

# 5. 按染色体位置排序，100-gene moving average
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gn <- as.data.frame(GenomicFeatures::genes(txdb))
# TxDb 的 gene ID 是 Entrez；表达矩阵行名是 HGNC symbol，必须先映射
sym <- AnnotationDbi::select(org.Hs.eg.db, keys = gn$gene_id,
                             keytype = "ENTREZID", columns = "SYMBOL")
pos <- data.frame(
  gene  = sym$SYMBOL,
  chr   = gn$seqnames[match(sym$ENTREZID, gn$gene_id)],
  start = gn$start[match(sym$ENTREZID, gn$gene_id)]
) %>%
  filter(!is.na(gene), gene %in% rownames(Er)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(chr, start)
stopifnot(nrow(pos) > 1000)
Er <- Er[pos$gene, , drop = FALSE]

moving_average <- function(mat, window = 100) {
  out <- mat
  for (ch in unique(pos$chr)) {
    idx <- which(pos$chr == ch)
    if (length(idx) == 0) next
    sub <- mat[idx, , drop = FALSE]
    for (i in seq_along(idx)) {
      lo <- max(1, i - floor(window / 2))
      hi <- min(length(idx), i + floor(window / 2))
      out[idx[i], ] <- colMeans(sub[lo:hi, , drop = FALSE])
    }
  }
  out
}
CNV_i <- moving_average(Er, window = 100)

# 6. 每个细胞中心化
CNV_c <- CNV_i - rep(colMeans(CNV_i), each = nrow(CNV_i))

# 7. baseline：hE02-hE06 非恶性上皮
baseline_clusters <- c("hE02_Enterocyte-GUCA2B", "hE03_Goblet-SPINK4",
                       "hE04_Goblet-REG4", "hE05_Stemlike-LGR5", "hE06_UnIdent")
base_idx <- which(sce$cluster_id %in% baseline_clusters)
if (length(base_idx) < 10) stop("baseline 上皮细胞不足；检查注释")
base_mat <- CNV_c[, base_idx, drop = FALSE]
BaseMax <- apply(base_mat, 1, max)
BaseMin <- apply(base_mat, 1, min)

CNV_f <- CNV_c
for (j in seq_len(ncol(CNV_f))) {
  x <- CNV_c[, j]
  y <- ifelse(x > BaseMax + 0.2, x - BaseMax,
              ifelse(x < BaseMin - 0.2, x - BaseMin, 0))
  CNV_f[, j] <- y
}
CNV_f[CNV_f > 1] <- 1; CNV_f[CNV_f < -1] <- -1
CNV_f[abs(CNV_f) < 0.1] <- 0

# 8. CNV score / correlation（按患者分别计算恶性参考；原文 Methods：same sample）
cnv_score <- colMeans(CNV_f^2)
stopifnot("patient_id" %in% colnames(colData(sce)))
cnv_cor <- rep(NA_real_, ncol(CNV_f))
is_malig <- rep(FALSE, ncol(CNV_f))
for (p in unique(as.character(sce$patient_id))) {
  idx_p <- which(as.character(sce$patient_id) == p)
  # 该患者的恶性候选：hC01–hC04；若不可用，用该患者 CNV score 前 5%
  malig_p <- idx_p[sce$cluster_id[idx_p] %in% c("hC01", "hC02", "hC03", "hC04")]
  if (length(malig_p) < 10) {
    malig_p <- idx_p[cnv_score[idx_p] >= quantile(cnv_score[idx_p], 0.95)]
  }
  if (length(malig_p) < 3) next
  avg_p <- rowMeans(CNV_f[, malig_p, drop = FALSE])
  cnv_cor[idx_p] <- apply(CNV_f[, idx_p, drop = FALSE], 2, function(x) cor(x, avg_p))
  is_malig[idx_p] <- cnv_score[idx_p] > quantile(cnv_score[idx_p], 0.9) &
    cnv_cor[idx_p] > 0.5
}
if (any(is.na(cnv_cor))) warning("部分患者 CNV correlation 为 NA；检查该患者恶性候选细胞数")

sce$cnv_score <- cnv_score
sce$cnv_correlation <- cnv_cor
sce$is_malignant <- is_malig
assay(sce, "cnv") <- CNV_f
saveRDS(sce, file.path(proj, "processed", "human_cnv.rds"))

write_tsv(data.frame(cell_id = colnames(sce), cluster_id = sce$cluster_id,
                     cnv_score = cnv_score, cnv_correlation = cnv_cor,
                     is_malignant = sce$is_malignant),
          file.path(proj, "results", "03_cnv_scores.tsv"))
```

> ⚠️ **导师批注**：CNV 的 baseline 依赖 `hE02–hE06` 的正确注释；如果这些簇在 Smart-seq2 中不存在或细胞太少，CNV 会失真。恶性判定不能只用 `CNV score`（高分也可能是增殖/高表达上皮），要结合 `CNV correlation`、上皮 marker（`EPCAM/KRT8/18`）和聚类。`infercnv` 可作敏感性分析，但不要把两套结果混在一起报。

### 4.6 注释与簇数核对

```r
#!/usr/bin/env Rscript
# code/03_06_validation.R
suppressPackageStartupMessages({ library(Seurat); library(readr); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
obj <- readRDS(file.path(proj, "processed", "human_annot.rds"))
meta <- obj@meta.data

counts <- meta %>% count(major_type, cluster_id, sort = TRUE)
write_tsv(counts, file.path(proj, "results", "03_cluster_counts.tsv"))

expected <- c(Myeloid = 13, ILC = 4, T = 18, B = 5)
observed <- table(meta$major_type)
print(rbind(expected = expected, observed = observed[names(expected)]))

# marker 方向检查
logtpm <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
check_marker <- function(cluster, genes, direction = "up") {
  g <- intersect(genes, rownames(logtpm))
  if (length(g) == 0) return(NA)
  m <- colMeans(logtpm[g, meta$cluster_id == cluster, drop = FALSE])
  other <- colMeans(logtpm[g, meta$cluster_id != cluster, drop = FALSE])
  if (direction == "up") mean(m > other) else mean(m < other)
}
stopifnot(
  check_marker("hM12_TAM-C1QC", c("C1QA","C1QB","C1QC","TREM2","MERTK","CD80")) > 0.5,
  check_marker("hM13_TAM-SPP1", c("SPP1","MARCO","VEGFA")) > 0.5
)
cat("[OK] 注释与 marker 方向核对通过\n")
```

---

## ⑤ 参数解释与调参空间

| 参数 | 原文取值 | 可选范围 | 改变后果 | 本项目裁决 |
|---|---|---|---|---|
| HVG 数量 | 约 2,000 | 1,000–3,000 | 少 → 稀有簇丢失；多 → 噪声 | 2,000（原文） |
| HVG 方法 | mean/dispersion | vst/mean.var.plot/dispersion | 影响簇结构 | `mean.var.plot`（最接近原文） |
| PCs | 10–25 | 8–30 | 少 → 欠分割；多 → 噪声 | 20（范围内），敏感性 10/25 |
| resolution | 0.6–2.0 | 0.4–2.5 | 直接改变簇数 | 1.0 主分析；扫描 0.6–2.0 |
| 聚类算法 | Louvain | Leiden/Louvain | 簇边界略有差异 | Louvain（Seurat 默认） |
| 第一轮 resolution | 文献未报告 | 0.4–0.8 | 影响大类纯度 | 0.6 |
| 相似性 alpha | 0.99 | 0.5–1.0 | 影响特征选择 | 0.99（原文） |
| CV folds | 10 | 5–10 | 影响 lambda 选择 | 10（原文） |
| 下采样最小簇 | 50 | 30–100 | 影响类别平衡 | 50（原文） |
| CNV window | 100 genes | 50–200 | 窗口小 → 噪声；大 → 平滑过度 | 100（原文） |
| CNV 表达过滤 | log2(mean TPM+1) < 3.5 | 3.0–4.0 | 影响基因数/CNV 分辨率 | 3.5（原文） |
| CNV baseline | hE02–hE06 | 所有非恶性上皮 | 影响 CNVf 零点 | 原文；敏感性用全部上皮 |
| 恶性判定 | CNV score + correlation | 单独 score / 单独 correlation | 影响恶性细胞数 | 两者结合 + 上皮 marker |

---

## ⑥ 预期结果与 QC 检查点

| 编号 | 指标 | 文献值 | 允许范围 | Go/No-Go | 我方值（S6 填） |
|---|---|---|---|---|---|
| Q09 | 10x 免疫簇数 | 38 | 37–39 | ✅/⚠️/❌ | |
| Q10 | Smart-seq2 免疫簇数 | 36 | 35–37 | ✅/⚠️/❌ | |
| Q11 | 人类髓系簇数 | 13 | 12–14 | ✅/⚠️/❌ | |
| Q12 | ILC/T/B 簇数 | 4/18/5 | ±1 | ✅/⚠️/❌ | |
| Q13 | 非免疫簇数 | 12 | 11–13 | ✅/⚠️/❌ | |
| Q14 | 恶性细胞判定 | 8 例肿瘤的恶性簇 | 方向一致；CNV score/correlation 分离 | ✅/⚠️ | |
| Q15 | hM12/hM13 marker 方向 | C1QC↑/SPP1↑ | 方向一致 | ✅/❌ | |

**Go/No-Go**：

- ✅ 13 个髓系簇全部存在（允许 ±1，但缺少 `hM09`/`hM04` 等关键簇要排查）。
- ✅ `hM12_TAM-C1QC` 的 `C1QA/B/C`、`TREM2`、`MERTK`、`CD80` 平均表达高于其他簇；`hM13_TAM-SPP1` 的 `SPP1`、`MARCO`、`VEGFA` 高于其他簇。
- ✅ 跨平台相似性矩阵中，同一标签的 10x/SS2 簇应对角/近对角；若大面积错位，检查基因 ID/标签。
- ✅ CNV 能把恶性与非恶性分开（`CNV score` 高 + `CNV correlation` 高），且 8 例肿瘤各自的恶性簇可识别。
- ❌ 若髓系簇 <11 或 >15：先查 QC/分辨率/注释，不要直接改标签。

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 簇数远多于 38/36 | resolution 过高；HVG 过多；批次未校正 | 降 resolution；检查 02 号 BF；用 marker 稳定性选参数 |
| 簇数远少于 38/36 | resolution 过低；HVG 过少；稀有簇被 QC 过滤 | 提高 resolution；检查 10x gene 下限；检查 DC marker |
| 髓系与 T 细胞混在一起 | 第一轮大类未分开；PCA 被淋巴主导 | 先粗分大类；必要时对髓系子集单独聚类 |
| `hM12`/`hM13` 分不开 | `C1QC`/`SPP1` dropout；10x 3′ 偏好 | 用 Smart-seq2 做 marker/DEG；10x 用于比例；检查 HVG 是否含这两个基因 |
| 自动注释标签错 | 最高分 marker 受批次/表达量影响 | 人工核对 marker 热图；用 `03_cluster_labels_final.tsv` 固定 |
| 跨平台相似性矩阵方向反 | 训练/测试写反 | 明确行=测试、列=训练；图注写清 |
| `cv.glmnet` 报 offset 长度错 | offset 未按样本重复 | `rep(log(f/(1-f)), n)`；检查 `n` |
| Ro/e 全为 1 | tissue 字段全是 Tumor | 从 metadata/文件名解析组织；检查 `table(cluster, tissue)` |
| CNV 全为 0 | baseline 与细胞同一分布；窗口过大 | 检查 baseline 簇；缩小窗口；检查 `CNV_c` 中心化 |
| CNV 把上皮都判恶性 | 只用 score，未用 correlation/上皮 marker | 两者结合；检查 `hE02–hE06` 是否被误判 |
| `TxDb` 基因名与矩阵不一致 | Ensembl vs Symbol | 统一 ID；用 `AnnotationDbi::select`；记录丢失 |
| 恶性簇患者特异 | 肿瘤异质性/CNV 差异 | 原文也是患者特异簇（hC01–hC04）；按患者分别判定，不要强行合并 |

---

## ⑧ 与文献图表的对应关系表

| 原文图/表 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 1A | `processed/human_smartseq2_clustered.rds`、`human_annot.rds` | Smart-seq2 热图 |
| Fig. 1B | `results/03_cluster_counts.tsv` | 两平台簇数比较 |
| Fig. 1C | `figures/`（CCA 可视化，可选） | 双平台联合 t-SNE；**只用于可视化** |
| Fig. S2A–S2D | `processed/human_*_clustered.rds` | 免疫亚群 t-SNE/热图 |
| Fig. S2E | `results/03_cross_platform_similarity.tsv`、`FigS2E_similarity.pdf/.png` | 跨平台相似性 |
| Fig. S2F–S2G | `results/03_cluster_markers.tsv` | B/ILC marker |
| Fig. S3A–S3C | `results/ro_e.tsv`、`Fig2C_ro_e.pdf/.png` | 平台/组织分布、Ro/e |
| Fig. S3D–S3F | `processed/human_cnv.rds`、`results/03_cnv_scores.tsv` | 非免疫 12 簇、CNV、恶性判定 |
| Fig. S3G–S3H | `results/03_cluster_markers.tsv` | 非免疫 marker、组织分布 |
| Fig. 2C | `results/ro_e.tsv` | 髓系组织富集 |
| Fig. S4A–S4B | `processed/human_annot.rds`（髓系子集） | 13 髓系簇 t-SNE/marker |

---

## 附录 A：最终标签表模板（`results/03_cluster_labels_final.tsv`）

```text
cluster_id	final_label	major_type	evidence	note
0	hM12_TAM-C1QC	Myeloid	C1QA/B/C, TREM2, MERTK, CD80	Fig. 2B/2D
1	hM13_TAM-SPP1	Myeloid	SPP1, MARCO, VEGFA	Fig. 2B/2D
...
```

> 落盘时不得保留 `<...>`；每一行必须有 marker 证据或跨平台相似性证据。

## 附录 B：簇数核对表（S6 填）

| 类别 | 文献值 | 10x 实测 | Smart-seq2 实测 | 差异原因 |
|---|---|---|---|---|
| 免疫总簇 | 38 / 36 | | | |
| 髓系 | 13 | | | |
| ILC | 4 | | | |
| T | 18 | | | |
| B | 5 | | | |
| 非免疫 | 12 | | | |
| 恶性（患者特异） | 4（hC01–hC04） | | | |

## 附录 C：CNV 参数速查

| 步骤 | 公式/参数 | 原文值 |
|---|---|---|
| 低表达过滤 | `log2(mean(TPM)+1) < 3.5` | 3.5 |
| 表达变换 | `E = log2(TPM/10 + 1)` | 10 |
| 基因中心化 | `Er = E − mean_i(E)` | — |
| 截断 | `[-3, 3]` | 3 |
| 滑动窗口 | 100 genes/chromosome | 100 |
| 细胞中心化 | `CNV_c = CNV_i − mean_j(CNV_i)` | — |
| baseline | hE02–hE06 | 5 簇 |
| 阈值 | `±0.2` | 0.2 |
| 画图 | `[-1,1]`；`[-0.1,0.1]→0` | 1 / 0.1 |
| 恶性分数 | `mean(CNV_f²)`；Pearson correlation | — |
