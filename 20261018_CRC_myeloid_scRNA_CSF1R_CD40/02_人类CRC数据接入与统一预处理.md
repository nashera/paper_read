# 02 人类 CRC 数据接入与统一预处理

> **目标文献**：Zhang L, Li Z, Skrzypczynska KM, et al. *Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer.* Cell. 2020;181(2):442–459.e29. DOI: [10.1016/j.cell.2020.03.048](https://doi.org/10.1016/j.cell.2020.03.048)
> **模块对应原文**：Results p.444（细胞数与平台比较）；STAR Methods「Single-cell RNA-Seq data processing」（p.e7–e8）；「Evaluation of batch effect with entropy metric」（p.e8）；「Batch effect correction and normalization」（p.e8）；「Cell doublet detection and removal」（p.e8）；Fig. S1A–S1F
> **前置依赖**：`01_环境搭建与数据获取.md`（`raw/human_crc_10x/`、`raw/human_crc_smartseq2/`、`sample_metadata.tsv`）；`PROJECT_SPEC.md` §4 契约（`processed/human_10x_qc.rds` 等）、§5 术语、§6 版本、§8 D1–D6
> **复现难度 / 预计工时 / 计算资源**：★★★☆☆ / 2–3 人日 / 32 GB 内存、8 核 CPU；Scrublet 与熵计算可并行

---

## ① 科学问题与逻辑位置

### 1.1 本模块要解决什么

把 `raw/` 里两种平台、两种定量口径、两种组织来源的矩阵，变成**同一套可比较的细胞级对象**：

- 统一基因 ID（HGNC symbol）与细胞元数据（患者、组织、平台、QC 指标）；
- 按原文阈值完成细胞过滤（Smart-seq2 的 median−3 MAD；10x 的 gene/UMI 区间；线粒体 10%/70%）；
- 用熵指标评估患者批次；
- 用 pagoda2 batch factor 做跨患者校正；
- 用 Scrublet + cluster-level 规则去除双细胞（尤其 `CD14⁺IL32⁺` 单核–NK 双细胞）；
- 输出 `processed/human_10x_qc.rds`、`human_smartseq2_qc.rds`、`human_10x_norm.rds`、`human_smartseq2_norm.rds`。

上游：01 号 `raw/**`、`sample_metadata.tsv`。
下游：`03_无监督聚类注释与恶性细胞CNV.md` 读取四个 `processed/*.rds`；04–07 号读取 `*_norm.rds`；14 号用 QC 指标核对 Q01–Q08。

### 1.2 两条平台链路的差异（必须先想清楚）

| 维度 | Smart-seq2 | 10x Chromium |
|---|---|---|
| 捕获 | 96 孔板、单细胞分选 | 微流控液滴、barcode |
| 定量 | 全长 cDNA → kallisto → tximport gene counts | Cell Ranger UMI counts |
| 长度偏好 | 无（全长） | 3′ 端 |
| 基因数 | 15,179 | 13,538 |
| 细胞数 | **10,468（已含 5,126 个既往 T 细胞）** | 43,817 |
| 主要用途 | DEG、GSVA、signature、LR | 聚类、稀有簇、比例 |
| 关键坑 | ERCC、GAPDH QC、无 UMI | 空滴、双细胞、dropout |

> **工程纪律**：两平台**分别 QC、分别归一化、分别聚类**；跨平台只在 03 号用 logistic regression 做"相似性映射"，**不合并矩阵做 DE**。原文也是这么做的。

---

## ② 背景原理

### 2.1 kallisto + tximport：Smart-seq2 的定量

- **kallisto**：伪比对（pseudoalignment）把 reads 映射到 transcriptome，用 EM 估计 transcript 丰度；`--bootstrap-samples=100` 给不确定性。
- **tximport**：把 transcript-level estimated counts 汇总到 gene-level；同时给出 `length`、`counts`、`abundance`（TPM）。
- **TPM**：
  $$
  \mathrm{TPM}_{ij}=\frac{10^6 \cdot \mathrm{RC}_{ij}/L_i}{\sum_{i'}\mathrm{RC}_{i'j}/L_{i'}}
  $$
  其中 `RC` 是 read count，`L_i` 是基因有效长度。TPM 先按长度归一，再按样本总量归一，**跨样本可比**。
- **log2(TPM+1)**：稳定方差，压缩高表达；原文后续分析用这个尺度。

### 2.2 Cell Ranger：10x 的定量

- `cellranger count`：比对（STAR）、过滤、UMI 校正、barcode 白名单、聚类（Graph-based）。
- 输出 `filtered_feature_bc_matrix`（barcode × gene UMI）。
- **UMI 不是 read count**：一个 mRNA 分子只计一次，减少 PCR 重复；但 3′ 端捕获、长度偏好、dropout 更明显。
- 原文用 Cell Ranger 2.1.0；若用新版，barcode 命名、`features.tsv` 列名、参考基因组版本都可能不同，必须记录。

### 2.3 QC 的四个维度

| 维度 | 指标 | 过低 | 过高 |
|---|---|---|---|
| **文库大小** | library size（总 counts/UMI） | 空滴、RNA 降解 | 双细胞、PCR 偏倚 |
| **基因数** | nGene（counts>0） | 低质量/空滴 | 双细胞、多核 |
| **线粒体** | percent_mito | — | 濒死/透化细胞 |
| **双细胞** | doublet score | — | 两个细胞一个液滴 |

原文的 Smart-seq2 规则：
- `library size < median − 3×MAD` 或 `nGene < median − 3×MAD` → 剔除；
- `percent_mito > 10%`（CD45⁺）或 `>70%`（CD45⁻）→ 剔除。

10x 规则：
- `nGene ∈ [600, 4000]`、`nUMI ∈ [1600, 25000]`；
- `percent_mito > 10%`（CD45⁺）或 `>70%`（CD45⁻）→ 剔除。

**MAD**：
$$
\mathrm{MAD}=\mathrm{median}\left(\left|x_i-\mathrm{median}(x)\right|\right)
$$
稳健标准差 ≈ `1.4826 × MAD`；`median − 3×MAD` 比固定阈值更适应患者/平台差异。

### 2.4 批次熵（entropy）

原文（Azizi et al. 2018）：
1. 取标准差最大的 1,000 个基因；
2. `Rtsne` 降到 2D；
3. 在 t-SNE 坐标上建 kNN 图（`k=30`，欧氏距离，`dbscan::kNN`）；
4. 对细胞 j，`p_j^t` = 30 邻域中来自患者 t 的比例；
5. `H_j = −Σ_t p_j^t log2 p_j^t`。

**直觉**：H 高 = 邻域混合均匀 = 批次弱；H 低 = 患者聚集 = 批次强。原文用 H 的分布比较 10x 与 Smart-seq2（Fig. S1D）。

### 2.5 pagoda2 batch factor

原文：
$$
\mathrm{RC}^{\text{nobatch}}_{ij}=\frac{C_{ij}}{\mathrm{BF}_{ik}}
$$
其中 `C_ij` 是基因 i 在细胞 j 的原始 count/UMI，`BF_ik` 是基因 i 在批次 k 的 batch factor。

**可审计的近似实现**（本项目主路径，若 pagoda2 可用则以其输出为准）：

1. 对每个基因 i、每个批次 k，计算批内均值 `μ_ik`；
2. 计算全局均值 `μ_i`；
3. `BF_ik = (μ_ik + ε) / (μ_i + ε)`；
4. `RC_nobatch = C / BF`。

> ⚠️ **导师批注**：pagoda2 内部对 BF 的估计包含方差/过度离散的处理，不是简单的均值比。上面的公式是**透明、可复现的近似**，用于在没有 pagoda2 的环境下跑通链路；正式结果应与 pagoda2 输出对比（至少在一个数据集上），差异写入 `logs/DECISIONS.md`（D5）。若用 pagoda2，代码中保留 `p2 <- pagoda2::Pagoda2$new(...)` 的入口并记录版本/commit；具体 API 需按当前版本核对。

### 2.6 TPM-like（10x）

10x UMI 没有长度信息；原文说"create a TPM-like value"。本项目实现：
$$
\mathrm{TPM\text{-}like}_{ij}=\frac{10^6 \cdot C_{ij}/L_i}{\sum_{i'} C_{i'j}/L_{i'}}
$$
`L_i` 从 `TxDb.Hsapiens.UCSC.hg38.knownGene` 或作者提供的基因长度表取得；没有长度时退化为 CPM（并记录）。最后 `log2(TPM-like + 1)`。

### 2.7 Scrublet 与 cluster-level 双细胞去除

- **Scrublet**：随机配对细胞表达向量，模拟双细胞；用 kNN 分类器给每个细胞打分；阈值由双峰分布自动确定。
- 原文还用 **cluster-level** 规则：删除 `CD14⁺IL32⁺` 簇——同时表达单核（`CD14`、`S100A8/9/12`）与细胞毒（`IL32`、`NKG7`、`GZMA/B`）的簇，典型单核–NK 双细胞。
- 复现时保留被删簇的细胞数、marker 均值、doublet score 均值，写入 `logs/02_doublet_removed.tsv`。

---

## ③ 原文逐句拆解

| 原文句（节选） | 中文转述 | 我的操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "Low-quality reads were filtered out when: (1) 'N' bases accounting for 10% read length; or (2) bases with quality < 5 account for 50% read length; or (3) containing adaptor sequences." | Smart-seq2 read 质控 | 若从 FASTQ 重跑才执行；主路径用 processed | 给了规则 | 不重跑；记录 |
| "The remaining paired-end reads were aligned to the UCSC hg38 ... using kallisto (version 0.43.1) with parameters '--bootstrap-samples=100 --threads=2'." | kallisto 参数 | 仅记录；processed matrix 已定量 | 给了版本/参数 | D1 |
| "tximport (version 1.9.12) was used to summarize ... gene-level counts." | tximport 汇总 | 若重跑则用；主路径检查作者矩阵 | 给了版本 | D1 |
| "Low-quality cells were discarded if the cell library size or the number of expressed genes ... smaller than ... medians of all cells minus 3 × median absolute deviation." | Smart-seq2 median−3 MAD | `02_02_qc_smartseq2.R` | 给了规则 | 照做 |
| "Cells were also removed if their proportions of mitochondrial gene expression were larger than 10% (70% for CD45− cells)." | 线粒体阈值 | 按 `CD45` 状态分支 | 给了阈值 | 照做；70% 仅 CD45⁻ |
| "Cell Ranger toolkit (version 2.1.0) ... GRCh38 and GRCm38." | 10x 上游 | 记录；processed matrix 已定量 | 给了版本 | D2 |
| "kept high quality cells with the detection threshold of 600–4,000 genes and 1,600–25,000 UMIs." | 10x 阈值 | `02_01_qc_human_10x.R` | 给了阈值 | 照做 |
| "Genes with the top 1000 standard deviation were used to map the data to a two-dimensional t-SNE by R package Rtsne (version 0.15)." | 熵的输入 | `compute_entropy()` | 给了 top 1000 + 版本 | 照做 |
| "k-NN graph (k = 30) ... using R function kNN from dbscan package (version 1.1.2)" | 熵的 kNN | `dbscan::kNN` | 给了 k=30 | 照做 |
| "H_j = −Σ_t p_j^t log2 p_j^t" | Shannon 熵 | 逐细胞计算 | 给了公式 | 照做 |
| "batch factor (BF) using the approach in pagoda2 ... RC_nobatch = C_ij / BF_ik" | batch factor | `estimate_batch_factor()`；pagoda2 可选 | 给了公式，未给 BF 估计细节 | D5：主路径用透明近似；与 pagoda2 对比 |
| "In Smart-seq2 data, TPM = ...; in 10x data, ... TPM-like ... log2(TPM+1)" | 归一化 | `normalize_counts()` | 给了公式 | 照做 |
| "We applied Scrublet ... default parameters ... removed the CD14+IL32+ cluster" | 双细胞 | `02_05_scrublet.py` | 给了默认参数 | 照做 |
| "5,126 single T cells ... combined" | 既往 T 细胞 | `02_03_merge_prior_T.R` | 给了数量 | 若文件缺失，记录 gap |

---

## ④ 逐步可运行代码

### 4.1 10x QC

```r
#!/usr/bin/env Rscript
# code/02_01_qc_human_10x.R
# 输入: raw/human_crc_10x/（10x mtx 目录，或作者提供的 matrix）
# 输出: processed/human_10x_qc.rds, results/02_10x_qc_metrics.tsv

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(Matrix)
  library(readr); library(dplyr); library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

# ---- 1. 读取矩阵（兼容 Read10X 目录与单个 mtx） ----
read_matrix <- function(dir) {
  if (file.exists(file.path(dir, "matrix.mtx")) ||
      file.exists(file.path(dir, "matrix.mtx.gz"))) {
    return(Seurat::Read10X(dir))
  }
  mtx <- list.files(dir, pattern = "\\.mtx(\\.gz)?$", full.names = TRUE)
  if (length(mtx) == 1) {
    m <- Matrix::readMM(mtx)
    barcodes <- list.files(dir, pattern = "barcodes", full.names = TRUE)
    features <- list.files(dir, pattern = "features|genes", full.names = TRUE)
    feat <- read.delim(features[1], header = FALSE, stringsAsFactors = FALSE)
    # 10x features.tsv 通常两列：Ensembl + symbol；genes.tsv 同
    rownames(m) <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
    colnames(m) <- readLines(barcodes[1])
    return(m)
  }
  stop("无法识别 10x 矩阵格式；请检查 raw/human_crc_10x/ 的文件")
}

counts <- read_matrix(file.path(proj, "raw", "human_crc_10x"))
# 若作者矩阵是 gene symbol，保持；若是 Ensembl，转 symbol 并记录丢失
if (all(grepl("^ENSG", rownames(counts)))) {
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(counts),
                               columns = "SYMBOL", keytype = "ENSEMBL")
  map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
  counts <- counts[map$ENSEMBL, , drop = FALSE]
  rownames(counts) <- map$SYMBOL
  write_tsv(map, file.path(proj, "logs", "02_ensembl_to_symbol.tsv"))
}

# ---- 2. 构建 Seurat 对象并计算 QC ----
obj <- CreateSeuratObject(counts = counts, project = "CRC_10x",
                          min.cells = 0, min.features = 0)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# 原文阈值
obj$pass_qc <- with(obj@meta.data,
  nFeature_RNA >= 600 & nFeature_RNA <= 4000 &
  nCount_RNA >= 1600 & nCount_RNA <= 25000 &
  percent.mt <= 10)

qc <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(patient_id = sub("_.*$", "", cell_id),
         tissue = "Tumor", platform = "10x") %>%
  select(cell_id, patient_id, tissue, platform, nGene = nFeature_RNA,
         nUMI = nCount_RNA, percent_mito = percent.mt, pass_qc)
write_tsv(qc, file.path(proj, "results", "02_10x_qc_metrics.tsv"))

obj <- subset(obj, subset = pass_qc)
sce <- as.SingleCellExperiment(obj)
sce$cell_id <- colnames(sce)
sce$patient_id <- sub("_.*$", "", colnames(sce))
sce$tissue <- "Tumor"   # 若 metadata 有 Tumor/Normal/Blood，按 metadata 覆盖
sce$platform <- "10x"
sce$nGene <- obj$nFeature_RNA
sce$nUMI <- obj$nCount_RNA
sce$percent_mito <- obj$percent.mt
sce$pass_qc <- TRUE

saveRDS(sce, file.path(proj, "processed", "human_10x_qc.rds"))
cat("[OK] 10x QC:", ncol(sce), "cells,", nrow(sce), "genes\n")
```

> ⚠️ **导师批注**：GEO 的 10x 文件可能是**每个患者一个目录**，也可能是合并矩阵。`read_matrix()` 目前只处理单目录；若为多患者目录，先循环读取、加患者前缀，再 `cbind`。`tissue` 字段必须从 metadata/文件名解析，不能默认全部 Tumor（上面的默认值只是占位，运行前必须按实际 metadata 覆盖）。

### 4.2 Smart-seq2 QC

```r
#!/usr/bin/env Rscript
# code/02_02_qc_smartseq2.R
# 输入: raw/human_crc_smartseq2/（counts + TPM + metadata）
# 输出: processed/human_smartseq2_qc.rds, results/02_ss2_qc_metrics.tsv

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(Matrix); library(readr); library(dplyr); library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

read_tsv_matrix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (is.matrix(obj)) return(obj)
    if (inherits(obj, "SummarizedExperiment")) return(as.matrix(assay(obj)))
    stop("RDS 对象既不是 matrix 也不是 SummarizedExperiment：", path)
  }
  df <- read.delim(path, row.names = 1, check.names = FALSE)
  as.matrix(df)
}

counts_file <- list.files(file.path(proj, "raw", "human_crc_smartseq2"),
                          pattern = "count", full.names = TRUE)[1]
tpm_file <- list.files(file.path(proj, "raw", "human_crc_smartseq2"),
                       pattern = "tpm|TPM", full.names = TRUE)[1]
stopifnot(!is.na(counts_file))

counts <- read_tsv_matrix(counts_file)
if (all(grepl("^ENSG", rownames(counts)))) {
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(counts),
                               columns = "SYMBOL", keytype = "ENSEMBL")
  map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
  counts <- counts[map$ENSEMBL, , drop = FALSE]
  rownames(counts) <- map$SYMBOL
}

# TPM：优先用作者提供的 TPM；否则用 counts + 基因长度计算
if (!is.na(tpm_file)) {
  tpm <- read_tsv_matrix(tpm_file)
  tpm <- tpm[intersect(rownames(counts), rownames(tpm)), , drop = FALSE]
  counts <- counts[rownames(tpm), , drop = FALSE]
} else {
  # 需 TxDb；若未安装，报错而不是猜
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    stop("缺少 TPM 且未安装 TxDb；请安装或提供作者 TPM 矩阵")
  }
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  gn <- GenomicFeatures::genes(txdb)
  len <- setNames(width(gn), names(gn))
  common <- intersect(rownames(counts), names(len))
  counts <- counts[common, , drop = FALSE]
  len <- len[common]
  rate <- sweep(counts, 1, len, "/")
  tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6
}

# QC 指标
lib_size <- colSums(counts)
n_gene <- colSums(counts > 0)
mito <- colSums(counts[grepl("^MT-", rownames(counts)), , drop = FALSE]) /
  pmax(lib_size, 1) * 100

# median − 3 MAD
mad_cut <- function(x) median(x) - 3 * mad(x, constant = 1)
pass <- lib_size > mad_cut(lib_size) &
        n_gene > mad_cut(n_gene) &
        mito <= 10   # CD45- 细胞改为 70

qc <- data.frame(cell_id = colnames(counts), patient_id = sub("_.*$", "", colnames(counts)),
                 tissue = "Tumor", platform = "Smart-seq2",
                 nGene = n_gene, nUMI = lib_size, percent_mito = mito,
                 pass_qc = pass)
write_tsv(qc, file.path(proj, "results", "02_ss2_qc_metrics.tsv"))

sce <- SingleCellExperiment(
  assays = list(counts = counts[, pass, drop = FALSE],
                tpm = tpm[, pass, drop = FALSE]))
sce$cell_id <- colnames(sce)
sce$patient_id <- sub("_.*$", "", colnames(sce))
sce$tissue <- "Tumor"
sce$platform <- "Smart-seq2"
sce$nGene <- n_gene[pass]
sce$nUMI <- lib_size[pass]
sce$percent_mito <- mito[pass]
sce$pass_qc <- TRUE
sce$prior_T_cell <- FALSE

saveRDS(sce, file.path(proj, "processed", "human_smartseq2_qc.rds"))
cat("[OK] Smart-seq2 QC:", ncol(sce), "cells,", nrow(sce), "genes\n")
```

> ⚠️ **导师批注**：CD45⁻ 细胞的 70% 线粒体阈值需要 `CD45` 状态；如果 processed matrix 的 metadata 没有 CD45 信息，先从 FACS 门控/作者 metadata 补，不能对所有细胞统一用 10% 或 70%。上面的 `mito <= 10` 只是 CD45⁺ 默认值，运行前必须按细胞类型分支。
>
> ⚠️ **计数口径（重要）**：GSE146771 的 Smart-seq2 processed matrix 已有 **10,468** 个细胞，原文明确这 10,468 是**合并 5,126 个既往 T 细胞之后**的总数（全文 3492–3494 与 3542–3543 行）。**不要**再额外加 5,126；4.3 只在核对发现矩阵确实不含既往 T 细胞时才执行。

### 4.3（可选）仅在矩阵不含既往 T 细胞时合并 5,126 个既往 T 细胞

```r
#!/usr/bin/env Rscript
# code/02_03_merge_prior_T.R
# 输入: processed/human_smartseq2_qc.rds + raw/human_crc_smartseq2/prior_t_cells_counts.tsv.gz
# 输出: processed/human_smartseq2_qc.rds（覆盖，含 prior_T_cell 标记）

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

sce <- readRDS(file.path(proj, "processed", "human_smartseq2_qc.rds"))

# 主路径：原文表格的 10,468 已含 5,126 个既往 T 细胞（全文 3492–3543 行）。
# 若当前对象已达 10,468 量级或已有 prior_T_cell 标记，直接退出，避免重复计数。
if (ncol(sce) >= 10000 || sum(sce$prior_T_cell, na.rm = TRUE) >= 5000) {
  message("Smart-seq2 对象已含既往 T 细胞（n=", ncol(sce), "），跳过合并")
  quit(save = "no")
}

prior_file <- file.path(proj, "raw", "human_crc_smartseq2",
                        "prior_t_cells_counts.tsv.gz")
if (!file.exists(prior_file)) {
  stop("缺少既往 T 细胞矩阵：", prior_file,
       "\n请从 Zhang et al. 2018 的补充数据获取；若无法获取，记录 gap 并跳过")
}
prior <- as.matrix(read.delim(prior_file, row.names = 1, check.names = FALSE))
common <- intersect(rownames(prior), rownames(sce))
prior <- prior[common, , drop = FALSE]
sce <- sce[common, , drop = FALSE]

# 仅保留 T 细胞（作者 metadata 标记；若无，则按 CD3D/CD3E/CD2 表达粗筛）
is_t <- rep(TRUE, ncol(prior))
if (file.exists(file.path(proj, "raw", "human_crc_smartseq2",
                          "prior_t_cells_metadata.tsv"))) {
  meta <- read_tsv(file.path(proj, "raw", "human_crc_smartseq2",
                             "prior_t_cells_metadata.tsv"))
  if ("cell_type" %in% colnames(meta)) {
    is_t <- meta$cell_type == "T"
  } else {
    warning("prior metadata 缺少 cell_type 列；保留全部 prior 细胞")
  }
}
prior <- prior[, is_t, drop = FALSE]

# 合并 counts（两套矩阵基因取交集；TPM 仅对原 sce 保留，prior 细胞 TPM 用 counts 重算）
merged <- cbind(assay(sce, "counts"), prior)
colnames(merged) <- make.unique(c(colnames(sce), colnames(prior)))
sce2 <- SingleCellExperiment(assays = list(counts = merged))
sce2$cell_id <- colnames(sce2)
sce2$patient_id <- c(sce$patient_id, rep("prior_Zhang2018", ncol(prior)))
sce2$tissue <- c(sce$tissue, rep("Tumor", ncol(prior)))
sce2$platform <- "Smart-seq2"
sce2$nGene <- colSums(merged > 0)
sce2$nUMI <- colSums(merged)
sce2$percent_mito <- colSums(merged[grepl("^MT-", rownames(merged)), , drop = FALSE]) /
  pmax(colSums(merged), 1) * 100
sce2$pass_qc <- TRUE
sce2$prior_T_cell <- c(rep(FALSE, ncol(sce)), rep(TRUE, ncol(prior)))

saveRDS(sce2, file.path(proj, "processed", "human_smartseq2_qc.rds"))
cat("[OK] 合并后 Smart-seq2:", ncol(sce2), "cells（既往 T:", ncol(prior), "）\n")
```

### 4.4 熵、batch factor、TPM/log2

```r
#!/usr/bin/env Rscript
# code/02_04_batch_norm.R
# 输入: processed/human_*_qc.rds
# 输出: processed/human_*_norm.rds, results/02_entropy_before.tsv, results/02_entropy_after.tsv

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(Rtsne); library(dbscan); library(Matrix); library(readr); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
set.seed(20261018)

# ---- 熵 ----
compute_entropy <- function(logmat, patient, k = 30, n_pcs = 2, seed = 20261018) {
  set.seed(seed)
  v <- apply(logmat, 1, sd)
  top <- order(v, decreasing = TRUE)[seq_len(min(1000, length(v)))]
  ts <- Rtsne::Rtsne(t(logmat[top, , drop = FALSE]), dims = n_pcs,
                     perplexity = min(30, floor((ncol(logmat) - 1) / 3)),
                     check_duplicates = FALSE)
  knn <- dbscan::kNN(ts$Y, k = k)
  # knn$id: cells x k
  H <- apply(knn$id, 1, function(idx) {
    p <- table(patient[idx]) / length(idx)
    -sum(p * log2(p))
  })
  H
}

# ---- batch factor（可审计的均值比近似） ----
estimate_batch_factor <- function(counts, batch, eps = 1e-6) {
  # counts: gene x cell; batch: factor
  genes <- rownames(counts)
  batches <- levels(factor(batch))
  mu_global <- rowMeans(counts) + eps
  BF <- matrix(1, nrow = nrow(counts), ncol = length(batches),
               dimnames = list(genes, batches))
  for (b in batches) {
    idx <- which(batch == b)
    mu_b <- rowMeans(counts[, idx, drop = FALSE]) + eps
    BF[, b] <- mu_b / mu_global
  }
  BF
}

# ---- TPM-like ----
tpm_like <- function(counts, gene_length = NULL) {
  if (is.null(gene_length)) {
    # 没有长度时退化为 CPM，并在日志中声明
    return(sweep(counts, 2, pmax(colSums(counts), 1), "/") * 1e6)
  }
  rate <- sweep(counts, 1, gene_length, "/")
  sweep(rate, 2, pmax(colSums(rate), 1), "/") * 1e6
}

process_one <- function(sce_path, out_path, tag) {
  sce <- readRDS(sce_path)
  counts <- as.matrix(assay(sce, "counts"))
  patient <- factor(sce$patient_id)
  logmat <- log2(tpm_like(counts) + 1)

  h_before <- compute_entropy(logmat, patient)
  write_tsv(data.frame(cell_id = colnames(sce), entropy = h_before),
            file.path(proj, "results", paste0("02_entropy_before_", tag, ".tsv")))

  BF <- estimate_batch_factor(counts, patient)
  counts_nb <- counts
  for (b in levels(patient)) {
    idx <- which(patient == b)
    counts_nb[, idx] <- counts[, idx, drop = FALSE] / BF[, b]
  }
  logmat_nb <- log2(tpm_like(counts_nb) + 1)
  h_after <- compute_entropy(logmat_nb, patient)
  write_tsv(data.frame(cell_id = colnames(sce), entropy = h_after),
            file.path(proj, "results", paste0("02_entropy_after_", tag, ".tsv")))

  sce2 <- sce
  assay(sce2, "logtpm") <- logmat_nb
  # BF 是 gene × batch，不能直接塞进 per-cell colData；
  # 这里存每细胞所属 batch 的 BF 中位数（摘要），完整矩阵另存。
  bf_med <- apply(BF, 2, median)
  colData(sce2)$batch_factor_median <- bf_med[as.character(patient)]
  saveRDS(BF, file.path(proj, "processed", paste0("02_batch_factor_", tag, ".rds")))
  saveRDS(sce2, out_path)
  cat("[OK]", tag, "熵中位数 before/after:",
      median(h_before), "/", median(h_after), "\n")
}

process_one(file.path(proj, "processed", "human_10x_qc.rds"),
            file.path(proj, "processed", "human_10x_norm.rds"), "10x")
process_one(file.path(proj, "processed", "human_smartseq2_qc.rds"),
            file.path(proj, "processed", "human_smartseq2_norm.rds"), "ss2")
```

> ⚠️ **导师批注**：`compute_entropy()` 对 43,817 个细胞做 `Rtsne` 可能很慢/占内存；可先按患者下采样（每患者 ≤2,000 细胞）计算熵，但要在 `logs/DECISIONS.md` 记录抽样方案。`estimate_batch_factor()` 是近似实现；若安装 pagoda2，请用其 batch factor 替换，并在 `results/` 保存两套 BF 的相关系数。

### 4.5 Scrublet 与 CD14⁺IL32⁺ 双细胞簇

```python
#!/usr/bin/env python3
# code/02_05_scrublet.py
"""Scrublet + cluster-level doublet removal for human 10x data.

输入: raw/human_crc_10x/ （10x mtx）
输出: processed/human_10x_scrublet.h5ad, results/02_doublet_scores.tsv
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

PROJ = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(PROJ, "raw", "human_crc_10x")
OUT = os.path.join(PROJ, "processed")
os.makedirs(OUT, exist_ok=True)

# 1. 读矩阵（10x mtx 目录）
adata = sc.read_10x_mtx(RAW, var_names="gene_symbols", cache=True)
adata.var_names_make_unique()

# 2. Scrublet（原文默认参数）
scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
scores, predicted = scrub.scrub_doublets()
adata.obs["doublet_score"] = scores
adata.obs["predicted_doublet"] = predicted

# 3. 聚类（用于 cluster-level 规则）
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=30, random_state=20261018)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20, random_state=20261018)
sc.tl.leiden(adata, resolution=1.0, random_state=20261018, key_added="leiden")

# 4. 找 CD14+IL32+ 双细胞簇：单核 + 细胞毒 marker 同时高
mono = ["CD14", "S100A8", "S100A9", "S100A12"]
cyto = ["IL32", "NKG7", "GZMA", "GZMB"]
def mean_expr(genes):
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return pd.Series(0.0, index=adata.obs_names)
    return pd.Series(
        np.asarray(adata[:, present].X.mean(axis=1)).ravel(),
        index=adata.obs_names)

adata.obs["mono_score"] = mean_expr(mono)
adata.obs["cyto_score"] = mean_expr(cyto)
cluster_stats = adata.obs.groupby("leiden", observed=True).agg(
    mono=("mono_score", "mean"),
    cyto=("cyto_score", "mean"),
    doublet=("doublet_score", "mean"),
    n=("doublet_score", "size"))
# 规则：mono 与 cyto 均高于全局中位数，且 doublet 均值高于全局 75 分位
mono_med = cluster_stats["mono"].median()
cyto_med = cluster_stats["cyto"].median()
d75 = cluster_stats["doublet"].quantile(0.75)
bad = cluster_stats.index[
    (cluster_stats["mono"] > mono_med) &
    (cluster_stats["cyto"] > cyto_med) &
    (cluster_stats["doublet"] > d75)]
print("CD14+IL32+ candidate clusters:", list(bad))
cluster_stats.to_csv(os.path.join(PROJ, "results", "02_doublet_cluster_stats.tsv"),
                     sep="\t")

# 5. 删除：predicted doublet 或 bad cluster
keep = (~adata.obs["predicted_doublet"].astype(bool)) & (~adata.obs["leiden"].isin(bad))
adata = adata[keep].copy()
adata.write_h5ad(os.path.join(OUT, "human_10x_scrublet.h5ad"))
adata.obs[["doublet_score", "predicted_doublet", "leiden"]].to_csv(
    os.path.join(PROJ, "results", "02_doublet_scores.tsv"), sep="\t")
print(f"[OK] retained {adata.n_obs} cells after doublet removal")
```

> ⚠️ **导师批注**：`sc.read_10x_mtx` 要求目录内有 `matrix.mtx(.gz)`、`barcodes.tsv(.gz)`、`features.tsv(.gz)`/`genes.tsv`。GEO 文件若命名不同，先重命名/转换，不要改 `read_10x_mtx` 的假设。`bad` 规则是**可审计的复现**，不是原文代码；原文只说"cluster-level approach"，未给阈值。把 `cluster_stats` 保存，便于敏感性分析。

### 4.6 QC 报告

```r
#!/usr/bin/env Rscript
# code/02_06_qc_report.R
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")

qc10 <- read_tsv(file.path(proj, "results", "02_10x_qc_metrics.tsv"))
qcss <- read_tsv(file.path(proj, "results", "02_ss2_qc_metrics.tsv"))
stopifnot(nrow(qc10) > 0, nrow(qcss) > 0)

p1 <- ggplot(qc10, aes(nGene, nUMI, color = pass_qc)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "10x QC", x = "nGene", y = "nUMI") + theme_bw()
p2 <- ggplot(qc10, aes(percent_mito, fill = pass_qc)) +
  geom_histogram(bins = 60) + labs(title = "10x mito%") + theme_bw()
p3 <- ggplot(qcss, aes(nGene, nUMI, color = pass_qc)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Smart-seq2 QC", x = "nGene", y = "library size") + theme_bw()
p4 <- ggplot(qcss, aes(percent_mito, fill = pass_qc)) +
  geom_histogram(bins = 60) + labs(title = "Smart-seq2 mito%") + theme_bw()

ggsave(file.path(proj, "figures", "FigS1B_qc_metrics.pdf"),
       p1 + p2 + p3 + p4, width = 12, height = 8)
ggsave(file.path(proj, "figures", "FigS1B_qc_metrics.png"),
       p1 + p2 + p3 + p4, width = 12, height = 8, dpi = 300)

summary_tbl <- bind_rows(
  qc10 %>% summarise(dataset = "10x", n_total = n(), n_pass = sum(pass_qc),
                     median_gene = median(nGene), median_umi = median(nUMI)),
  qcss %>% summarise(dataset = "Smart-seq2", n_total = n(), n_pass = sum(pass_qc),
                     median_gene = median(nGene), median_umi = median(nUMI))
)
write_tsv(summary_tbl, file.path(proj, "results", "02_qc_summary.tsv"))
print(summary_tbl)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 原文取值 | 可选范围 | 改变后果 | 本项目裁决 |
|---|---|---|---|---|
| Smart-seq2 低质量阈值 | median − 3 MAD | 2–4 MAD | 更严 → 细胞少/稀有簇丢失；更松 → 低质量 | 3 MAD（原文） |
| Smart-seq2 mito | 10%（CD45⁺）/70%（CD45⁻） | 5–20%/50–80% | 影响免疫 vs 肿瘤细胞保留 | 原文；按 CD45 分支 |
| 10x nGene | 600–4,000 | 500–6,000 | 更严 → DC/稀有簇丢失 | 原文 |
| 10x nUMI | 1,600–25,000 | 1,000–30,000 | 更严 → 低深度细胞丢失 | 原文 |
| 10x mito | 10% | 5–20% | 更严 → 濒死细胞少 | 原文 |
| 熵 top 基因 | 1,000 SD | 500–2,000 | 影响熵的稳定性 | 1,000（原文） |
| 熵 k | 30 | 15–50 | 邻域大小影响 H 绝对值 | 30（原文） |
| Rtsne perplexity | 文献未报告 | 20–50 | 影响 2D 结构 | `min(30, (n-1)/3)`，敏感性 20/50 |
| BF 估计 | pagoda2 | 均值比/中位数比/ComBat | 影响归一化后的聚类 | D5：均值比近似；与 pagoda2 对比 |
| Scrublet rate | 默认 | 0.01–0.1 | 影响双细胞数 | 原文默认；0.06 为常见值（需按数据核对） |
| 双细胞簇规则 | 原文只说 cluster-level | 阈值可调 | 误删真实细胞/漏删双细胞 | 透明规则 + 保存 cluster stats |

---

## ⑥ 预期结果与 QC 检查点

| 编号 | 指标 | 文献值 | 允许范围 | Go/No-Go | 我方值（S6 填） |
|---|---|---|---|---|---|
| Q01 | 10x 细胞数 | 43,817 | 41,626–46,008 | ✅/⚠️/❌ | |
| Q02 | Smart-seq2 细胞数 | 10,468 | 9,945–10,991 | ✅/⚠️/❌ | |
| Q03 | 10x 基因数 | 13,538 | 12,861–14,215 | ✅/⚠️/❌ | |
| Q04 | Smart-seq2 基因数 | 15,179 | 14,420–15,938 | ✅/⚠️/❌ | |
| Q05 | 既往 T 细胞 | 5,126 | 4,870–5,382 | ✅/⚠️ | |
| Q06 | 双细胞簇 `CD14⁺IL32⁺` | 去除 | 是 | ✅/❌ | |
| Q07 | 熵（每患者） | 文献未报告绝对值 | 与 Fig. S1D 同量级；after ≥ before | ✅/⚠️ | |
| Q08 | 归一化后无负值/NA | 文献未报告 | `summary(logtpm)` 无 NA，范围合理 | ✅/❌ | |

**Go/No-Go**：

- ✅ QC 后细胞数在 ±5%；双平台基因数与原文一致（±5%）。
- ✅ `CD14⁺IL32⁺` 簇被识别并删除；删除的细胞数与 doublet score 记录在案。
- ✅ 熵中位数在 BF 校正后不下降（允许 ±0.05 波动）；若下降，检查 BF 方向是否写反。
- ❌ 若 Smart-seq2 细胞数 < 9,000 或 10x < 40,000：先查阈值/文件，不要直接放宽阈值凑数（L0 §12）。

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 10x 读取失败 | GEO 文件名不是标准 10x 三件套 | 重命名/转换；或写自定义 reader；不要改 `Read10X` 的期望 |
| 所有细胞 mito=0 | 基因名是小写 `mt-` 或 Ensembl | 统一 symbol；人用 `^MT-`，鼠用 `^mt-` |
| 基因数骤降 | Ensembl→Symbol 映射丢失/重复 | 保存 `logs/02_ensembl_to_symbol.tsv`；检查 `length(intersect())` |
| Smart-seq2 细胞数远低于 10,468 | 矩阵是 TPM 而非 counts；MAD 计算在 TPM 上 | 确认 `counts` assay；MAD 只对 library size/nGene（counts）做 |
| 熵校正后更低 | BF 方向写反（乘 instead of 除） | `RC_nobatch = C / BF`；检查 `BF` 是否 batch/global |
| pagoda2 装不上 | GitHub 依赖 | 用本文件的均值比实现；记录 D5 |
| Scrublet 报 `expected_doublet_rate` 错误 | 10x 细胞数/双细胞率不匹配 | 设 0.06 或按 10x 手册估计；记录 |
| 双细胞簇规则删掉真实单核 | mono/cyto 阈值过松 | 检查 cluster stats；用 marker 表达和 doublet score 双重确认 |
| 合并既往 T 细胞后基因数变化 | 两套矩阵基因交集 | 取交集；记录丢失基因；不要 fill 0 强行并集 |
| `Rtsne` 内存爆 | 43,817 细胞 × 1,000 基因 | 按患者下采样；或用 `openTSNE`；记录抽样 |
| `logtpm` 出现 Inf | TPM 计算除零/长度为 0 | 加 `pmax(...,1)`；过滤长度 0 的基因 |

---

## ⑧ 与文献图表的对应关系表

| 原文图/表 | 本模块的对应产物 | 说明 |
|---|---|---|
| Fig. S1A | `logs/00_download_manifest.tsv`、`results/02_qc_summary.tsv` | 实验设计与细胞数 |
| Fig. S1B | `figures/FigS1B_qc_metrics.pdf/.png` | 基因数/核糖体比例（核糖体比例可选） |
| Fig. S1C | `results/02_qc_summary.tsv`（基因类别计数） | 细胞因子/CD/配体/受体/TF 检出数 |
| Fig. S1D | `results/02_entropy_before_*.tsv`、`02_entropy_after_*.tsv` | 批次熵直方图 |
| Fig. S1E | `results/02_dropout_ratio.tsv`（可选） | dropout 比值；需按基因计算 |
| Fig. S1F | 可选 KEGG 富集 | 10x 特异 dropout 基因 |
| Fig. 1A–1C | `processed/human_*_norm.rds` | 后续 03 号聚类输入 |
| STAR Methods 数据表 | `results/02_qc_summary.tsv` | 细胞数/基因数核对 |

---

## 附录 A：QC 指标定义与计算口径

| 指标 | 定义 | 尺度 | 备注 |
|---|---|---|---|
| `nGene` | 该细胞中 counts > 0 的基因数 | 个 | 10x 用 `nFeature_RNA` |
| `nUMI` | 该细胞的总 counts/UMI | 个 | Smart-seq2 是 read counts |
| `percent_mito` | 线粒体基因 counts / 总 counts × 100 | % | 人 `^MT-`，鼠 `^mt-` |
| `library size` | 总 counts | 个 | Smart-seq2 的 `nUMI` 实际是 library size |
| `entropy` | 30 邻域患者分布的 Shannon 熵 | bit | 越高批次越弱 |
| `BF` | 基因 × 批次的乘性因子 | 无量纲 | `RC_nobatch = C / BF` |
| `doublet_score` | Scrublet 的双细胞概率 | 0–1 | 阈值由双峰分布定 |

---

## 附录 B：输出契约核对

| 产物 | 本模块是否生成 | 下游 |
|---|---|---|
| `processed/human_10x_qc.rds` | ✅ `02_01` | 03, 04, 06, 07 |
| `processed/human_smartseq2_qc.rds` | ✅ `02_02` + `02_03` | 03, 04, 06, 07 |
| `processed/human_10x_norm.rds` | ✅ `02_04` | 03, 04, 05, 06, 07 |
| `processed/human_smartseq2_norm.rds` | ✅ `02_04` | 03, 04, 05, 06, 07 |
| `results/02_10x_qc_metrics.tsv` | ✅ | 14 |
| `results/02_ss2_qc_metrics.tsv` | ✅ | 14 |
| `results/02_qc_summary.tsv` | ✅ | 14 |
| `results/02_entropy_before_*.tsv` | ✅ | 14 |
| `results/02_entropy_after_*.tsv` | ✅ | 14 |
| `figures/FigS1B_qc_metrics.pdf/.png` | ✅ | 00, 14 |
