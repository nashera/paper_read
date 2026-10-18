# 07 PDX 体内模型与单细胞图谱

> **目标文献**：França GS, Baron M, King BR, Bossowski JP, Bjornberg A, Pour M, Rao A, Patel AS, Misirlioglu S, Barkley D, Tang KH, Dolgalev I, Liberman DA, Avital G, Kuperwaser F, Chiodin M, Levine DA, Papagiannakopoulos T, Marusyk A, Lionnet T, Yanai I. *Cellular adaptation to cancer therapy along a resistance continuum.* **Nature (2024) 631:876–890.** DOI: 10.1038/s41586-024-07690-9
> **模块对应原文**：Fig. 4a–g；Extended Data Fig. 9a–g；Methods「Generation of patient derived xenografts and drug treatments」「Tumour dissociation for scRNA-seq」「Tissue staining」「10X Genomics library preparation and sequencing」
> **前置依赖**：`03`/`04`（状态与模块）、`05`（inferCNV 流程）、`01`（GSE247690 数据）
> **复现难度 / 预计工时 / 计算资源**：★★★★☆（计算）/ 2–3 人日；体内实验另需动物伦理与数月周期

---

## ① 科学问题与逻辑位置

**问题**：

1. 体外定义的 persister 与 adapted 状态，是否在体内出现？
2. 剂量梯度（persistence）与长期治疗（resistance）产生的细胞状态是否不同？
3. 耐药肿瘤内部是否存在状态异质性（同一瘤内 SOX17^high 与 CD44^high 共存）？

**逻辑位置**：这是"从培养皿到动物"的桥接模块。它把 04 号的模块 A–F 作为**参照字典**，投射到 PDX 单细胞图谱上，检验体外结论的体内相关性。

**上游**：`raw/geo/GSE247690/*`（persistence 10 样本、resistance 7 样本）、`results/04_gene_modules.csv`、`results/05_infercnv.rds`（参考流程）。
**下游**：`results/07_pdx_clusters.csv`、`results/07_pdx_module_overlap.csv`、`figures/Fig4d–g_*`。

---

## ② 背景原理

### 2.1 PDX 模型

- 患者来源异种移植（PDX）保留原发肿瘤的组织结构与部分异质性，比细胞系更接近临床。
- 原文使用 **BRCA2 突变 HGSOC PDX（ID 66799）**，NSG 小鼠，皮下植入乳腺脂肪垫。
- **Persistence 实验**：固定治疗窗口内比较 4 个剂量（D1–D4），观察剂量依赖的转录状态。
- **Resistance 实验**：0.33 mg/kg 连续三轮；第一轮出现缓解→停药→复发→再治疗无反应（R1），提示耐药；R1* 传代形成 R2、R3。

### 2.2 为什么用 fastMNN

PDX 样本跨批次、跨小鼠、跨测序运行，批次效应强。`fastMNN`（batchelor）用互近邻（MNN）对齐，保留生物异质性；原文对 persistence 数据使用 fastMNN。Resistance 数据未明确说明，本项目对两套数据分别评估批次后决定（D9）。

### 2.3 物种混合去除

PDX 肿瘤含小鼠基质/免疫细胞。原文用 **CellRanger multi-species detection**，只保留 `human` 细胞。公开 metadata 已是 G1 人细胞；仍需核对 barcode 与 counts。

### 2.4 模块投射

把 04 号模块的基因作为签名，在 PDX 每个 cluster 上计算平均表达/模块评分，再用**超几何检验**（GeneOverlap）评估 cluster marker 与模块基因的重叠。原文报告 cluster 9 与模块 C/D/E 显著重叠。

### 2.5 免疫荧光

SOX17（谱系因子）与 CD44（间质/干性）双染，DAPI 核染。预期：V1（载体）以 SOX17^high CD44^low 为主；R3.2 出现 SOX17^low CD44^high 区域，与 scRNA-seq 的异质性一致。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 操作 | 原文参数 | 复现决策 |
|---|---|---|---|---|
| "BRCA2-mutant HGSOC patient-derived xenograft (PDX) model" | BRCA2 突变 HGSOC PDX | 模型 | ID 66799 | 分析公开数据 |
| "tumours were randomized ... 100–200 mm³" | 肿瘤 100–200 mm³ 开始 | 分组 | — | 方案记录 |
| "D1 0.13, D2 0.23, D3 0.33, D4 0.66 mg/kg" | 四个剂量 | 给药 | 5 天/周，灌胃 | 固定 |
| "collected after 50–57 days of treatment" | 50–57 天取瘤 | 时间点 | — | 记录 |
| "humane endpoint (2,000 mm³)" | 人道终点 | 伦理 | 2,000 mm³ | 固定 |
| "three passages ... 0.33 mg/kg" | 三轮耐药 | 传代 | — | 方案记录 |
| "drug holiday in the first passage" | 第一轮停药期 | 设计 | — | 固定 |
| "scRNA-seq ... only cells assigned to 'human'" | 只保留人细胞 | 物种过滤 | — | 公开 metadata 已过滤 |
| "fastMNN method for correction" | 批次校正 | fastMNN | — | 主口径（D9） |
| "cluster 9 ... overlap with modules C, D and E" | cluster 9 映射 | 超几何 | — | 以重叠而非编号判定 |
| "CD44 and SOX17 ... immunofluorescence" | 双染验证 | IF | 10 μg/ml | 方案记录 |
| "extensive CNAs relative to persister" | cluster 9 CNA 多 | inferCNV | — | 复用 05 |

**原文未报告**：PDX 的传代数/代数稳定性；小鼠性别/年龄（Methods 有 NSG 6 周雌性）；取瘤的精确时间点（50–57 天为区间）；IF 的定量方法；cluster 数的选择依据。均登记 D 表。

---

## ④ 逐步可运行代码

### 4.1 读取 persistence / resistance 数据

```r
# code/07_01_pdx.R
suppressPackageStartupMessages({
  library(Seurat); library(batchelor); library(SingleCellExperiment)
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
})
set.seed(20261018)

read_pdx <- function(counts_path, meta_path) {
  meta <- read.delim(meta_path, check.names = FALSE)
  stopifnot("cell" %in% colnames(meta))
  rownames(meta) <- make.unique(meta$cell)
  if ("orig.ident" %in% colnames(meta) && !"sample" %in% colnames(meta)) {
    meta$sample <- meta$orig.ident
  }
  if ("percent.mito" %in% colnames(meta)) {
    meta$percent.mito <- as.numeric(meta$percent.mito) * 100
  }
  mat <- as.matrix(read.csv(counts_path, row.names = 1, check.names = FALSE))
  if (mean(colnames(mat) %in% rownames(meta)) < 0.5) mat <- t(mat)
  common <- intersect(colnames(mat), rownames(meta))
  stopifnot(length(common) / nrow(meta) > 0.95)
  list(mat = mat[, common, drop = FALSE], meta = meta[common, , drop = FALSE])
}

per <- read_pdx(
  "raw/geo/GSE247690/GSE247690_pdx_persistence_counts_all.csv.gz",
  "raw/geo/GSE247690/GSE247690_pdx_persistence_metadata_g1.tsv.gz")
res <- read_pdx(
  "raw/geo/GSE247690/GSE247690_pdx_resistance_counts_all.csv.gz",
  "raw/geo/GSE247690/GSE247690_pdx_resistance_metadata_g1.tsv.gz")
```

### 4.2 建对象、fastMNN、聚类

```r
build_pdx <- function(x, use_mnn = TRUE) {
  srt <- CreateSeuratObject(counts = x$mat, meta.data = x$meta)
  srt <- NormalizeData(srt, scale.factor = 1e4, verbose = FALSE)
  srt <- FindVariableFeatures(srt, nfeatures = 2000, verbose = FALSE)
  srt <- ScaleData(srt, verbose = FALSE)
  srt <- RunPCA(srt, npcs = 30, verbose = FALSE)

  if (use_mnn) {
    # batchelor::fastMNN 需要 matrix/SCE，不能直接传 Seurat 对象
    objs <- lapply(SplitObject(srt, split.by = "sample"),
                   function(x) as.SingleCellExperiment(x))
    mnn <- do.call(fastMNN, c(objs, list(d = 15, k = 20, approximate = TRUE)))
    srt[["mnn"]] <- CreateDimReducObject(
      embeddings = reducedDim(mnn, "corrected"),
      key = "MNN_", assay = DefaultAssay(srt))
    srt <- FindNeighbors(srt, reduction = "mnn", dims = 1:15, verbose = FALSE)
  } else {
    srt <- FindNeighbors(srt, dims = 1:10, verbose = FALSE)
  }
  srt <- FindClusters(srt, resolution = 0.5, algorithm = 1, verbose = FALSE)
  srt <- RunUMAP(srt, reduction = if (use_mnn) "mnn" else "pca",
                 dims = 1:15, verbose = FALSE)
  srt
}

per_srt <- build_pdx(per, use_mnn = TRUE)   # 原文明确 fastMNN
res_srt <- build_pdx(res, use_mnn = TRUE)   # D9：同样校正，另做未校正敏感性

# 契约产物：cell, experiment, sample, cluster
bind_rows(
  data.frame(cell = colnames(per_srt), experiment = "persistence",
             sample = per_srt$sample, cluster = as.character(per_srt$seurat_clusters)),
  data.frame(cell = colnames(res_srt), experiment = "resistance",
             sample = res_srt$sample, cluster = as.character(res_srt$seurat_clusters))
) %>% write.csv("results/07_pdx_clusters.csv", row.names = FALSE)
```

### 4.3 聚类组成与 UMAP

```r
plot_pdx <- function(srt, tag) {
  p1 <- DimPlot(srt, group.by = "sample", reduction = "umap") +
    ggtitle(paste0(tag, " — by sample"))
  p2 <- DimPlot(srt, group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0(tag, " — clusters"))
  ggsave(paste0("figures/", tag, "_umap_sample.pdf"), p1, width = 6, height = 5)
  ggsave(paste0("figures/", tag, "_umap_clusters.pdf"), p2, width = 6, height = 5)

  comp <- as.data.frame(table(srt$sample, srt$seurat_clusters))
  colnames(comp) <- c("sample", "cluster", "n")
  comp <- comp %>% group_by(sample) %>%
    mutate(fraction = n / sum(n)) %>% ungroup()
  write.csv(comp, paste0("results/07_", tag, "_cluster_composition.csv"),
            row.names = FALSE)
  comp
}
comp_per <- plot_pdx(per_srt, "07_pdx_persistence")
comp_res <- plot_pdx(res_srt, "07_pdx_resistance")
```

### 4.4 cluster marker 与模块重叠

```r
# 每个数据集分别找 marker
markers_per <- FindAllMarkers(per_srt, test.use = "MAST",
  logfc.threshold = 0.25, min.pct = 0.2, min.cells.feature = 10,
  only.pos = TRUE, verbose = FALSE)
markers_res <- FindAllMarkers(res_srt, test.use = "MAST",
  logfc.threshold = 0.25, min.pct = 0.2, min.cells.feature = 10,
  only.pos = TRUE, verbose = FALSE)
write.csv(markers_per, "results/07_pdx_persistence_markers.csv", row.names = FALSE)
write.csv(markers_res, "results/07_pdx_resistance_markers.csv", row.names = FALSE)

# 与模块 A–F 超几何重叠
library(GeneOverlap)
mod <- read.csv("results/04_gene_modules.csv")
overlap_res <- lapply(split(mod$gene, mod$module), function(gs) {
  lapply(split(markers_res$gene, markers_res$cluster), function(mg) {
    go <- newGeneOverlap(mg, gs, genome.size = 20000)
    go <- testGeneOverlap(go)
    data.frame(p = getPval(go), odds = getOddsRatio(go),
               n_overlap = length(intersect(mg, gs)))
  }) %>% bind_rows(.id = "cluster")
}) %>% bind_rows(.id = "module")
write.csv(overlap_res, "results/07_pdx_module_overlap.csv", row.names = FALSE)

# 以重叠而非 cluster 编号判定"体内适应状态"
best <- overlap_res %>% filter(p < 0.05) %>% arrange(cluster, p)
print(best)
```

### 4.5 关键基因热图（Fig. 4f）

```r
genes <- c("SOX17","PAX8","WT1","KRT8","CDKN2A","IFI6","IFI27","IFITM3",
           "VIM","CD44","SMAD3","DDIT4","CITED2","SGK1","PGK1","PFKP",
           "ACLY","SCD","MBOAT7","SREBF1","CTPS1","TYMS","GPX4","STEAP3",
           "FTL","SOD2","GGH","GSTM3","PRDX1","CYP1B1","TXNRD1","NQO1",
           "SQSTM1","NRF2")
genes <- intersect(genes, rownames(res_srt))
library(ComplexHeatmap)
avg <- AverageExpression(res_srt, features = genes,
                         group.by = "seurat_clusters", slot = "data")$RNA
pdf("figures/Fig4f_pdx_markers.pdf", width = 8, height = 10)
Heatmap(t(scale(t(as.matrix(avg)))), cluster_rows = TRUE, cluster_columns = TRUE,
        name = "scaled mean", column_title = "PDX cluster")
dev.off()
```

### 4.6 体内 CNA（Extended Data Fig. 9g）

```r
# 复用 05 号模块的 inferCNV 流程；参考 = 各自 vehicle (V1/V2)
# 注意：PDX 的参考应为 vehicle 肿瘤，而不是细胞系 C
# 运行后比较 cluster 9 与其他 cluster 的 CNA 负荷
cna_burden <- colSums(abs(den@expr.data - 1) > 0.1)
boxplot(cna_burden ~ res_srt$seurat_clusters, las = 2,
        ylab = "CNA burden (fraction of genes altered)",
        main = "PDX resistance: CNA burden by cluster")
```

### 4.7 免疫荧光定量（方案）

```text
切片：5–7 μm；H&E 相邻切片
固定：1.6% PFA；透化：3% BSA + 0.1% Triton X-100，1 h
一抗：SOX17（R&D AF1924-SP，10 μg/ml）、CD44（R&D MAB7045-SP，10 μg/ml），3 h RT
二抗：Alexa 594/488（1:200）；DAPI 1 μg/ml
成像：Keyence BZ-X；每样本 3 张切片
定量：Fiji/QuPath → 每视野 SOX17+、CD44+ 细胞比例；DAPI 分割 → 共表达
统计：按视野（n≥5）做混合效应模型或 Wilcoxon；报告 mean ± SD
```

---

## ⑤ 参数解释与调参空间

| 参数 | 原文 | 本项目 | 改变后果 |
|---|---|---|---|
| 动物 | NSG，6 周雌性 | 方案记录 | 免疫缺陷程度影响植入 |
| 植入细胞数 | 1×10^6 + Matrigel 1:1 | 固定 | 过少成瘤率低 |
| 给药 | 5 天/周，灌胃 | 固定 | 暴露曲线不同 |
| 肿瘤起始体积 | 100–200 mm³ | 固定 | 过小变异大 |
| 人道终点 | 2,000 mm³ | 固定 | 伦理要求 |
| 聚类 resolution | 未报告 | 0.5 | 影响 cluster 数 |
| 批次校正 | fastMNN（persistence） | 主口径 | 不校正会按样本分离 |
| 物种过滤 | human only | 公开 metadata | 需核对 |
| 模块重叠检验 | 超几何 | GeneOverlap | 基因背景大小影响 P |
| IF 定量 | 未报告 | QuPath/Fiji | 需固定阈值并盲法 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 原文/预期 | 判据 |
|---|---|---|
| 样本数 | persistence 10；resistance 7 | 与 GEO 一致 |
| persistence UMAP | 剂量组混合 | 无按剂量的强分离 |
| resistance UMAP | 出现耐药特异 cluster | R2/R3 富集 |
| 未处理样状态 | cluster 0/6 类似体外 C | 模块 A 高 |
| persister 样 cluster | cluster 4 类似 P10/P320 | 模块 A 高 + 部分 E |
| 适应样 cluster | cluster 9 类似模块 C/D/E | 超几何 P<0.05 |
| 关键基因 | CD44/SMAD3/DDIT4/PFKP/CTPS1/TYMS/GPX4↑；SOX17/PAX8/WT1/CDKN2A↓ | 方向一致 |
| CNA | cluster 9 CNA 负荷高 | 与 persister 比较 |
| IF | R3.2 异质性 | SOX17^low CD44^high 区域存在 |
| 批间 | 两个 vehicle 可整合 | 无显著批次分离 |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| UMAP 按样本分离 | 批次效应 | fastMNN；检查 `SplitObject` 的样本名 |
| 小鼠基因混入 | 物种过滤失败 | 只保留 human metadata；检查基因名前缀 |
| cluster 编号与原文不同 | 聚类随机性/参数 | 用模块重叠而非编号报告 |
| 重叠检验全不显著 | marker 基因太少/背景大小 | 用 top markers；genome.size 固定 20,000 |
| CDKN2A 检测不到 | 表达低/ dropout | 用模块评分或免疫组化补充 |
| IF 背景高 | 一抗浓度/透化不足 | 滴定抗体；加 blocking |
| 成瘤率低 | 细胞活力/Matrigel 比例 | 新鲜细胞；1:1 比例 |
| 肿瘤生长曲线差异大 | 植入位置/小鼠批次 | 随机分组；同批小鼠 |
| 取瘤时间不一致 | 人道终点 | 记录并作为协变量 |

---

## ⑧ 与文献图表的对应关系

| 原文图 | 内容 | 本项目产物 | 说明 |
|---|---|---|---|
| Fig. 4a | persistence 实验设计 | `figures/Fig4a_persistence_scheme.pdf/png` | 方案图 |
| Fig. 4b | persistence UMAP + 组成 | `figures/07_pdx_persistence_umap_*.pdf/png` | fastMNN |
| Fig. 4c | resistance 实验设计 | `figures/Fig4c_resistance_scheme.pdf/png` | 三轮 |
| Fig. 4d | resistance UMAP + 组成 | `figures/07_pdx_resistance_umap_*.pdf/png` | |
| Fig. 4e | cluster × 模块重叠 | `results/07_pdx_module_overlap.csv` | 超几何 |
| Fig. 4f | 关键基因热图 | `figures/Fig4f_pdx_markers.pdf/png` | 34 基因 |
| Fig. 4g | IF（SOX17/CD44） | `figures/Fig4g_if_*.pdf/png` | 方案 + 定量 |
| Extended Data Fig. 9a–d | 生长曲线 | `results/07_pdx_growth_curves.csv` | 方案/数字化 |
| Extended Data Fig. 9e | H&E + IF 全片 | `figures/FigED9e_if_panel.pdf/png` | |
| Extended Data Fig. 9f | 剂量依赖基因热图 | `figures/FigED9f_dose_genes.pdf/png` | Mfuzz 可用于分簇 |
| Extended Data Fig. 9g | PDX CNA | `figures/FigED9g_pdx_cna.pdf/png` | inferCNV |

> **交付判定**：两个实验的 cluster 可解释、至少一个 cluster 与模块 C/D/E 显著重叠、关键基因方向一致、IF 方案完整；不要求 cluster 编号与原文一致。
