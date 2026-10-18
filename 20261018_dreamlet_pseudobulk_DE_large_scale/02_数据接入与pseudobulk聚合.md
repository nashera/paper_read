# 02 数据接入与 pseudobulk 聚合

> **目标文献**：Hoffman GE, et al. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* Nature Communications (2026) 17:9597. DOI: 10.1038/s41467-026-75680-8
> **模块对应原文**：Fig. 1（workflow 第一步）；Fig. 2A–B（pseudobulk 的内存与 CPU 时间）；Methods「Efficient computation of pseudobulk using on-disk memory」「Performance comparison for computing pseudobulk」「Creating pseudobulk from a large dataset…」
> **前置依赖**：`01_环境搭建与数据获取.md`（`raw/` 三个 h5ad、`logs/01_sessionInfo.txt`）
> **复现难度 / 预计工时 / 计算资源**：★★★☆☆ / 1–2 人日 / 内存 16–32 GB（on-disk 后端）；若用 in-memory 稀疏矩阵需 ≥64 GB

---

## ① 科学问题与逻辑位置

- **上游**：`raw/*.h5ad`（细胞 × 基因原始计数 + 细胞元数据）。
- **下游**：`processed/pb_<dataset>.rds`——**list**：`$pb`（`SingleCellExperiment`，每个细胞类型一个 `genes × samples` 的计数 assay）、`$cellCounts`（**样本 × 细胞类型**的细胞数矩阵，03 号第一级权重的输入）、`$config`；同时落盘 `processed/id_map_<dataset>.tsv`（Ensembl→SYMBOL，D2）。
- **本模块回答的问题**：**"把单细胞计数压成样本级 pseudobulk"这一步，怎么做才既不爆内存、又不丢信息、还能在大数据上快一个数量级？**

这一步是整条流水线的"地基"：如果 pseudobulk 的样本/细胞类型对错了，后面所有 DE、VP、GSA 都会安静地错（L0 附录 C「静默失败」）。

---

## ② 背景原理

### 2.1 pseudobulk 的数学定义

设细胞 $c$ 属于样本 $s$（`sample_id`）与细胞类型 $k$（`cluster_id`），基因 $g$ 的原始计数为 $y_{c,g}$。定义

$$
Y_{s,k,g} \;=\; \sum_{c \in \mathcal{C}_{s,k}} y_{c,g},
\qquad
L_{s,k} \;=\; \sum_{g} Y_{s,k,g},
\qquad
\mathcal{C}_{s,k} = \{c:\ \text{sample}(c)=s,\ \text{cluster}(c)=k\}.
$$

- $Y$ 是**计数之和**，仍是非负整数；**不做跨样本归一化**（归一化在 03 号）。
- $L_{s,k}$ 是该 pseudobulk 的文库大小。
- 若某样本在某细胞类型中没有细胞，则该格为缺失（不是 0）；dreamlet 用 `min.cells` 决定是否保留该样本。

**为什么求和而不是求平均？** 因为下游模型是"负二项/泊松式的计数生成过程 + log 变换"。求和保留了"观测到多少 reads"这一信息，使 03 号的 Poisson 权重有意义；求平均会丢失文库大小，且对低表达基因产生大量非整数噪声。

### 2.2 内存的量级：为什么必须 on-disk

原文给了一组关键算术（Methods）：

- 18,000 基因 × 2,000,000 细胞 × 8 字节（double）= **288 GB**（dense）。
- 若 80% 为零，`sparseMatrix` 只存非零项：约 $0.2 \times 18{,}000 \times 2{,}000{,}000 \times 8$ 字节 ≈ **58 GB**（仅数值，不含索引；若按 12 B/非零项计入索引则约 86 GB——原文的 58 GB 用的是 8 B 口径）。
- 仍不可行 → 数据留在磁盘，用 **HDF5/H5AD + DelayedArray** 按块读取。

一般化：

$$
\text{dense} = G \times C \times 8\ \text{B},\qquad
\text{sparse} \approx (1-\text{sparsity}) \times G \times C \times 12\ \text{B}.
$$

> ⚠️ **导师批注**：这里的 12 B = 8 B 数值 + 4 B 索引，是常见 CSR/CSC 的粗估；实际还取决于压缩与索引宽度。**复现性能时要用实测峰值内存，不要用这个公式交差。**

### 2.3 DelayedArray / beachmat 的工作方式

```text
H5AD (HDF5)
   │  zellkonverter::readH5AD(use_hdf5 = TRUE)
   ▼
DelayedMatrix（只是一个"读取计划"，不占内存）
   │  aggregateToPseudoBulk()
   ▼
按块（block）读取 → 块内按 (cluster, sample) 求和 → 累加到 genes × samples 矩阵
   │  beachmat / DelayedMatrixStats / sparseMatrixStats
   ▼
pseudobulk（小矩阵，可入内存）
```

- **块大小**：`h5adBlockSizes = 1e9`（约 1 GB）是 `aggregateToPseudoBulk()` 的默认自动块大小。块越大越快、内存越高。
- **beachmat**：Bioconductor 的 C++ 矩阵访问 API，为 DelayedMatrix / sparseMatrix / HDF5Matrix 提供统一的、按块的 C++ 接口；这是"同样算法、快一个数量级"的工程根源。
- **muscat 的 `aggregateData()`** 用的是 `scuttle::summarizeAssayByGroup()`，对 in-memory 矩阵没问题，但对 on-disk 矩阵会反复触发小块读取 → 慢。

### 2.4 为什么"按样本 × 细胞类型"而不是"按样本"

跨受试者 DE 的核心是**细胞类型特异**：同一疾病在 microglia 与神经元中的效应不同。按 `(sample, cluster)` 聚合保留了这一层；若只按样本聚合，就把细胞类型混合了，得到的只是"组织整体"的效应。

### 2.5 pseudobulk 的代价（必须写进讨论）

聚合是**信息有损**的：

1. 丢失样本内的表达分布（双峰、burst、离群细胞）——原文 Discussion 明确承认。
2. 若某细胞类型在样本中只有 1–2 个细胞，pseudobulk 的抽样噪声极大；`min.cells` 是对策。
3. 细胞类型注释错误会被"平均"进 pseudobulk，且无法在样本级纠正。

这不是缺陷，而是**用"跨受试者推断"换"细胞级异质性"** 的自觉取舍。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "Dreamlet creates pseudobulk data from the raw read counts for each sample and cell cluster stored in a SingleCellExperiment object." | 从 SCE 的原始计数按 sample × cluster 聚合 | `aggregateToPseudoBulk(assay="counts", sample_id=, cluster_id=)` | 未给 `fun` | 默认 `fun="sum"`；D 表无争议 |
| "Storing 18 K genes across 2 M cells with double precision would require 288 Gb memory." | dense 内存估算 | 体检 | 给了数字 | — |
| "If the dataset is so sparse that 80% of the entries are zero, loading the entire dataset as a sparseMatrix object requires 58 Gb." | 稀疏也放不下 | 用 on-disk | 给了数字 | — |
| "zellkonverter … uses a DelayedArray backend to provide a seamless interface to an on-disk H5AD dataset." | zellkonverter + DelayedArray | `readH5AD(use_hdf5=TRUE)` | 给了包名 | — |
| "we use the beachmat library to summarize each block containing all genes and a subset of cells into pseudobulk." | 用 beachmat 分块求和 | `aggregateToPseudoBulk()` 内部 | 未给块大小 | 默认 1 GB；大内存可调到 2–4 GB |
| "For eight cell types and 1000 donors, this reads data into memory in 8000 chunks to keep memory usage low." | naive 法读 8000 块 | 08 号对照实现 | 给了 8000 | — |
| "aggregateToPseudoBulk() consolidates these into fewer operations." | dreamlet 合并读取 | 同上 | — | — |
| "Each method was run ten times … average memory usage and CPU time … 512 Gb memory … one thread." | 基准设置 | 08 号 | 给了 | D13 |

---

## ④ 逐步可运行代码

### 4.1 读入并清洗（以 PsychAD 为例）

```r
# code/02_00_read_sce.R
# 输入：raw/psychad_ad/PsychAD_r0_Dec_28_2022.h5ad
# 输出：processed/sce_psychad_raw.rds（可选的中间层；大文件可跳过，直接进 4.2）
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
})

h5 <- "raw/psychad_ad/PsychAD_r0_Dec_28_2022.h5ad"
sce <- readH5AD(h5, use_hdf5 = TRUE, verbose = TRUE)

# 统一 assay 名：zellkonverter 常把第一个 assay 命名为 "X"
assayNames(sce)[1] <- "counts"

# 关键断言：counts 必须是非负整数
.stopif <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
.stopif(all(assay(sce, "counts") >= 0), "counts 含负值")
.stopif(is(assay(sce, "counts"), "integer") ||
        all(assay(sce, "counts") == round(assay(sce, "counts"))),
        "counts 非整数：可能读到了归一化矩阵（见 01 号 4.4）")

# 分组因子（术语表口径）
sce$Dx <- factor(sce$Dx_AD, levels = c("Control", "AD"))
sce$subclass <- factor(sce$subclass)
sce$SubID <- factor(sce$SubID)
sce$Channel <- factor(sce$Channel)
sce$poolID <- factor(sce$poolID)

saveRDS(sce, "processed/sce_psychad_raw.rds")
```

### 4.2 聚合为 pseudobulk

```r
# code/02_01_pseudobulk.R
# 用法：Rscript code/02_01_pseudobulk.R psychad
# 输出：processed/pb_<dataset>.rds
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dreamlet)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[1] else "psychad"

cfg <- list(
  psychad = list(
    h5 = "processed/sce_psychad_raw.rds",
    cluster_id = "subclass",   # 22 类
    sample_id  = "Channel"     # 一次 10x 建库
  ),
  kfoury = list(
    h5 = "processed/sce_kfoury_raw.rds",
    cluster_id = "cells",
    sample_id  = "ID"
  ),
  covid = list(
    h5 = "processed/sce_covid_raw.rds",
    cluster_id = "major_subset",
    sample_id  = "scRNASeq_sample_ID"
  )
)[[dataset]]

sce <- readRDS(cfg$h5)

# 去掉关键元数据缺失的细胞（否则聚合会得到 NA 组）
keep <- !is.na(sce[[cfg$cluster_id]]) & !is.na(sce[[cfg$sample_id]])
if (sum(!keep) > 0) message("drop cells with NA id: ", sum(!keep))
sce <- sce[, keep]
sce[[cfg$cluster_id]] <- droplevels(factor(sce[[cfg$cluster_id]]))
sce[[cfg$sample_id]]  <- droplevels(factor(sce[[cfg$sample_id]]))

# 聚合（原始计数求和）
t0 <- Sys.time()
pb <- aggregateToPseudoBulk(
  sce,
  assay      = "counts",
  cluster_id = cfg$cluster_id,
  sample_id  = cfg$sample_id,
  fun        = "sum",
  scale      = FALSE,
  BPPARAM    = SerialParam(progressbar = TRUE),
  h5adBlockSizes = 1e9,
  verbose    = TRUE
)
message("pseudobulk elapsed: ", format(Sys.time() - t0))

# 质量断言：非空、整数、列名可解释
stopifnot(ncol(pb) > 0, nrow(pb) > 0)
stopifnot(all(assay(pb, 1) >= 0))
stopifnot(all(assay(pb, 1) == round(assay(pb, 1))))

# ---------------------------------------------------------------------------
# 基因 ID 统一（D2）：Ensembl → SYMBOL，唯一口径在这里落地
# 之后所有文档（05/06/07/11/12）只消费 SYMBOL
# ---------------------------------------------------------------------------
ids <- rownames(pb)
is_ensembl <- grepl("^ENSG", ids[1])

if (is_ensembl) {
  suppressPackageStartupMessages({ library(org.Hs.eg.db); library(scater) })
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = ids,
                               keytype = "ENSEMBL", columns = "SYMBOL")
  sym <- map$SYMBOL[match(ids, map$ENSEMBL)]
  mapped <- !is.na(sym)
  sym[!mapped] <- ids[!mapped]                 # 未映射的保留原 ID
  new_ids <- scater::uniquifyFeatureNames(ids, sym)  # 重复 symbol 追加 Ensembl 后缀
  stopifnot(!anyDuplicated(new_ids))
  rownames(pb) <- new_ids
} else {
  mapped <- rep(TRUE, length(ids))
  sym <- ids
}

id_map <- data.frame(original_id = ids, symbol = sym, mapped = mapped,
                     final_id = rownames(pb))
write.table(id_map, sprintf("processed/id_map_%s.tsv", dataset),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("ID map: ", sum(!mapped), " unmapped, ",
        sum(duplicated(sym[mapped])), " duplicated symbols")

# 关键基因存在性断言（防止后面 PTPRG/MT-/MALAT1 静默取空）
if (dataset == "psychad") {
  stopifnot("PTPRG" %in% rownames(pb), "PDE10A" %in% rownames(pb))
}
if (dataset == "covid") {
  stopifnot(any(grepl("^MT-", rownames(pb))), "MALAT1" %in% rownames(pb))
}

# 每样本每类型细胞数（03 号第一级权重的输入）
cc <- cellCounts(pb)
saveRDS(list(pb = pb, cellCounts = cc, config = cfg),
        sprintf("processed/pb_%s.rds", dataset))

# 数据字典
write.csv(as.data.frame(cc), sprintf("results/02_cellcounts_%s.csv", dataset))
```

### 4.3 逐数据集的最小检查（必须做）

```r
# code/02_02_check_pb.R
pb_obj <- readRDS("processed/pb_psychad.rds")
pb <- pb_obj$pb; cc <- pb_obj$cellCounts

cat("assays (cell types):", length(assayNames(pb)), "\n")
print(assayNames(pb))

# 每个细胞类型的 样本数 × 基因数
info <- data.frame(
  celltype = assayNames(pb),
  n_samples = vapply(assayNames(pb), function(k) ncol(assay(pb, k)), 1L),
  n_genes   = vapply(assayNames(pb), function(k) nrow(assay(pb, k)), 1L),
  median_cells = vapply(assayNames(pb), function(k) median(cc[, k]), 1),
  min_cells    = vapply(assayNames(pb), function(k) min(cc[, k]), 1)
)
print(info, row.names = FALSE)
write.csv(info, "results/02_pb_summary_psychad.csv", row.names = FALSE)

# 术语表口径的标签检查（防大小写/单复数错配）
stopifnot(setequal(assayNames(pb), c(
  "EN_L2_3_IT","EN_L3_5_IT_1","EN_L3_5_IT_2","EN_L3_5_IT_3","EN_L5_6_NP",
  "EN_L6_CT","EN_L6_IT","EN_NF","IN_ADARB2","IN_LAMP5","IN_PVALB","IN_PVALB_CHC",
  "IN_SST","IN_VIP","Oligo","OPC","Astro","Micro_PVM","CD8_T","PC","VLMC","Endo"
)))
```

### 4.4 性能实测（为 08 号提供数据）

```r
# code/02_03_bench_pseudobulk.R
# 只测 aggregateToPseudoBulk 的 wall time 与 peak memory（用 gc() 前后差 + /usr/bin/time 更准）
library(dreamlet); library(SingleCellExperiment); library(BiocParallel)

bench_once <- function(sce, cluster_id, sample_id, block = 1e9) {
  gc(reset = TRUE)
  t0 <- Sys.time()
  pb <- aggregateToPseudoBulk(sce, assay = "counts", cluster_id = cluster_id,
                              sample_id = sample_id, fun = "sum",
                              h5adBlockSizes = block, verbose = FALSE)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mem <- sum(gc()[, "max used"])   # 单位 MB；仅作粗估，精确值请用外部 /usr/bin/time -v
  list(pb = pb, seconds = dt, mem_mb = mem)
}
```

> 精确峰值内存请用 shell 的 `/usr/bin/time -v`（Linux）或 `gtime -v`（macOS）包住整条 `Rscript` 命令，而不是在 R 内部估。

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| `assay` | `"counts"` | 归一化 assay（禁止） | 用归一化值聚合 → 计数模型失效，p 值全错 |
| `fun` | `"sum"` | `"mean"`, `"median"`, `"num.detected"` | 非 `sum` 时不能用 dreamlet 的计数权重；`num.detected` 可做"检出率"分析但非本复现 |
| `scale` | `FALSE` | `TRUE`（CPM 缩放） | `TRUE` 会改变文库大小语义；本复现保持 FALSE |
| `cluster_id` | 各队列注释列 | 更粗/更细的注释 | 决定"细胞类型特异"的粒度；不能与原文 22 类混用 |
| `sample_id` | AD `Channel`；Kfoury `ID`；COVID `scRNASeq_sample_ID` | 受试者 ID（错误） | 用受试者 ID 会把技术重复合并 → 丢失重复测量结构，VP/随机效应失效 |
| `h5adBlockSizes` | `1e9`（1 GB） | `2e8`–`4e9` | 越大越快、峰值内存越高；>4 GB 在 32 GB 机器上危险 |
| `BPPARAM` | `SerialParam()`（基准）/ `SnowParam(8)`（生产） | 任意 | 并行加速聚合；内存近似线性增长 |
| `min.cells`（下游） | 见 03/05 | 3/5/10 | 过滤掉细胞数过少的样本；过严会丢样本、过松会引入噪声 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期 | 允许波动 | Go/No-Go |
|---|---|---|---|
| AD pseudobulk 列数（样本） | 586 | ±5% | ❌ 查 `Channel` 定义与 QC |
| AD 细胞类型数 | 22 | 完全一致 | ❌ 标签拼写 |
| AD microglia 的样本数 | **不设 586 的硬预期**（`min.cells=5` 会过滤掉细胞数少的样本）；记录实际值，与作者代码/补充图对照 | — | ⚠️ 记录 |
| COVID pseudobulk 列数 | 110（去重后每位 donor 一个样本） | 完全一致（D10） | ❌ 查样本选择 |
| Kfoury 列数 | 样本数（含 Benign/Distal/Involved/Tumor 四类） | — | ❌ 查 `Status` |
| 每个 assay 的值 | 非负整数 | 必须 | ❌ |
| `cellCounts` 的 NA | 0 | 必须 | ❌ |
| 聚合耗时（AD，on-disk） | 与 08 号基准同量级；本机 1–3 h | — | ⚠️ 只比相对值 |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 某些细胞类型 assay 缺失 | 该类型细胞数 < `min.cells` 或标签为空 | `table(sce[[cluster_id]], useNA="ifany")`；检查注释流程 |
| `aggregateToPseudoBulk` 报 "no method for DelayedMatrix" | `beachmat`/`DelayedMatrixStats` 版本不匹配 | `BiocManager::install(c("beachmat","DelayedMatrixStats","sparseMatrixStats"))` |
| 内存爆掉 | 用了 `use_hdf5=FALSE` 或 block 太大 | 改 `use_hdf5=TRUE`；`h5adBlockSizes=5e8` |
| pseudobulk 出现大量 0 列 | `sample_id` 与 `cluster_id` 组合不存在 | 交叉表检查；可能 `sample_id` 取错列（如取了 `donor_id`） |
| 数值与作者结果差 2 倍 | 作者用了 `scale=TRUE` 或 CPM | 对照作者脚本 `aggregateToPseudoBulk(...)` 的实参 |
| 聚合后基因名重复 | h5ad 的 `var_names` 有重复 SYMBOL | 用 Ensembl 作行名；或 `make.unique()` 并记录 |
| 并行后结果列顺序变化 | `SnowParam` 任务返回顺序 | 聚合结果按列名对齐；不要依赖列顺序 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 1（workflow 第一格） | `processed/pb_<dataset>.rds` | 每个细胞类型一个 assay |
| Fig. 2A（peak memory） | `results/benchmark_memory.csv`（08 号汇总） | 本模块提供单次测量函数 |
| Fig. 2B（CPU time） | `results/benchmark_time.csv` | 同上 |
| Fig. 2C 的输入 | 各数据集的 `pb_*.rds` | DE 基准的起点 |
| Fig. 5C / 3B / 4B 的均值–方差趋势 | `pb_*.rds` → 03 号 | 趋势图来自 03 号的 `plotVoom()` |

> **交付物**：`processed/pb_psychad.rds`、`pb_kfoury.rds`、`pb_covid.rds`；`results/02_pb_summary_*.csv`；`results/02_cellcounts_*.csv`。完成后进入 `03_两阶段精度权重与均值方差趋势.md`。
