# 06 ATAC-seq 染色质可及性与转录因子 motif

> **目标文献**：França GS, Baron M, King BR, Bossowski JP, Bjornberg A, Pour M, Rao A, Patel AS, Misirlioglu S, Barkley D, Tang KH, Dolgalev I, Liberman DA, Avital G, Kuperwaser F, Chiodin M, Levine DA, Papagiannakopoulos T, Marusyk A, Lionnet T, Yanai I. *Cellular adaptation to cancer therapy along a resistance continuum.* **Nature (2024) 631:876–890.** DOI: 10.1038/s41586-024-07690-9
> **模块对应原文**：Fig. 3b–f；Extended Data Fig. 8a–f；Methods「ATAC-seq sample preparation and sequencing」「ATAC-seq data processing and analysis」
> **前置依赖**：`01_数据与代码可得性、环境搭建.md`（GSE247685、参考基因组、HOMER）、`03`/`04`（分选亚群定义）、`05`（CNA 区段）
> **复现难度 / 预计工时 / 计算资源**：★★★★★ / 3–5 人日 / 64 GB RAM；若从 SRA 重跑需 ≥500 GB 磁盘与 8–16 核

---

## ① 科学问题与逻辑位置

**问题**：

1. 沿连续谱，染色质可及性是否**累积**变化，而非随机波动？
2. 哪些转录因子（TF）的 motif 可及性随适应逐步关闭/打开？
3. 撤药后染色质状态是否稳定？再暴露能否"恢复"？
4. 这些变化是否只是 CNA 的副产物？

**逻辑位置**：04 号模块给出转录程序，05 号模块给出遗传背景；本模块提供**表观遗传层**的证据，是"epigenetically reinforced cell states"与"cellular memory"的核心。

**上游**：`raw/geo/GSE247685/*`（**20 样本计数**的 consensus 表 + 6 个系列级 bigWig；28 样本的 fastq 在 SRA）、`results/05_cna_segments.csv`。
**下游**：`results/06_diff_peaks.csv`、`results/06_motif_enrichment.csv`、`figures/Fig3c–f_*`。

---

## ② 背景原理

### 2.1 ATAC-seq 与 OMNI-ATAC

ATAC-seq 用工程化 Tn5 转座酶插入开放染色质，测序读段两端即开放区域边界。**OMNI-ATAC** 通过去垢剂（digitonin、Tween-20）去除线粒体与背景，降低假信号，适合冷冻/分选样本。

### 2.2 peak calling → IDR → consensus

- **MACS2**：按读段堆积 call peak；`--nomodel` 常用于 ATAC。
- **IDR（Irreproducible Discovery Rate）**：两个生物学重复的 peak 排序一致性；IDR<0.05 保留。
- **consensus peaks**：合并所有样本的可重复 peak，形成统一坐标，便于跨样本计数。

### 2.3 差异可及性与归一化

DESeq2 用负二项模型处理 count 数据；ATAC 的 count 表由 deepTools 在 consensus peaks 上统计。归一化考虑文库深度；差异判据 adj P<0.001 且 FC≥3。

### 2.4 motif 富集

HOMER `findMotifsGenome.pl` 对差异 peak 区域做已知 motif 富集（超几何检验），背景为同长度随机基因组区域或指定背景。motif 可及性变化提示**上游 TF 活性**变化，但不是直接结合证据。

### 2.5 "记忆"的操作化

- T320：亲本适应系（在 320 μM 维持）。
- T320_nodrug：撤药 10 天（原文 "parental T320 grown without drug for 10 days"）。
- T320_drug：撤药后再暴露 320 μM。
- T320W：撤药 **120 天**。
- T320W_drug：T320W 再暴露 320 μM。

比较"更开放 peak"的信号：T320W 撤药时下降、再暴露时恢复 → 染色质记忆。

### 2.6 数据缺口与两条路线

| 路线 | 输入 | 能做什么 | 不能做什么 |
|---|---|---|---|
| **A. 公开 consensus counts（默认，严格）** | `GSE247685_consensus_peaks_kuramochi_all.tsv.gz`（chr/start/end + **20 样本计数**） | Fig. 3c–e、PCA、motif、CNA 控制 | **不含 T320/T320W 加药/撤药 8 样本**，无法做 Fig. 3f |
| **B. SRA 重跑（补充）** | BioProject PRJNA1040167 的 fastq | Fig. 3f 的 8 个样本；或从零重建 peak/IDR | 算力与磁盘成本高 |
| **C. 系列级 bigWig** | 6 个 `.bw` | 信号可视化（如 IGV/热图） | 不能替代计数矩阵做 DE |

**关键事实**（2026-10-18 实测）：`GSE247685_consensus_peaks_kuramochi_all.tsv.gz` 的第一行是
`chr,start,end,C_1,T5_1,T10_1,T20_1,T40_1,T80_1,T160_1,T320_1,C_2,T5_2,T10_2,T20_2,T40_2,T80_2,T160_2,T320_2,P10_1,P10_2,P320_1,P320_2`。
即 **20 个样本的计数已经公开**，可直接跑 DESeq2；缺的是 `T320_nodrug/T320_drug/T320W_nodrug/T320W_drug` 8 个样本。SOFT 中各样本 `supplementary_file_1 = NONE`，6 个 bigWig 是系列级轨道，不能按样本拆分。D7 据此更新。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 操作 | 原文参数 | 复现决策 |
|---|---|---|---|---|
| "modified version of the OMNI-ATAC protocol" | 改良 OMNI-ATAC | 建库 | — | 路线 A |
| "approximately 50,000 cells per condition ... nuclei extraction" | 5 万细胞 | 核提取 | 50,000 | 固定 |
| "transposition ... 37°C for 40 min ... 1,000 rpm" | 转座条件 | — | 40 min | 固定 |
| "TrimGalore 0.5.0 (-q 20, --length 30)" | 去接头 | 质控 | q20/30 | 固定 |
| "Bowtie2 v2.3.4.1 (--very-sensitive, -k 10, -X 2000)" | 比对 | — | 参数 | 固定 |
| "PicardTools 2.18.2 MarkDuplicates" | 去重复 | — | — | 固定 |
| "SAMTools 1.9 ... properly mapped pairs" | 过滤 | — | — | 固定 |
| "downsampled ... to match the lowest depth" | 降采样 | — | — | 固定 |
| "MACS2 v2.2.4 ... blacklisted regions were excluded" | call peak | — | blacklist | 固定 |
| "IDR 2.0.3 ... IDR score <0.05" | 重复一致性 | — | 0.05 | 固定 |
| "HOMER v4.15 ... consensus peak set" | 合并 peak | — | — | 固定 |
| "deepTools 3.1 ... read count table" | 计数 | — | — | 固定 |
| "DESeq2 1.18.1 ... adjusted P < 0.001 and fold change ≥ 3" | 差异 peak | — | 0.001/3 | 固定 |
| "top 2,000 most variable peaks ... PCA" | PCA | — | 2,000 | 固定 |
| "Zavit62 ... sort the mean peak signal" | 排序 | — | — | 固定 |
| "HOMER's findMotifsGenome ... known transcription factor motifs" | motif | — | known | 固定 |
| "permutation tests (1,000) ... regioneR" | CNA 重叠 | — | 1,000 | 固定 |
| "removing ... peaks that overlap with gained or deleted copy-number regions" | CNA 控制 | — | — | 固定 |
| "bamCoverage --binSize 1 --normalizeUsing RPGC" | 信号轨 | — | — | 固定 |
| "computeMatrix --binSize 1 --referencePoint center ... plotHeatmap" | 热图 | — | — | 固定 |
| "Differences in peak signal ... Wilcoxon tests" | 统计 | — | — | 固定 |

**原文未报告**：MACS2 的具体 `--nomodel/--shift/--extsize`；IDR 的输入排序方式；HOMER 的背景选择；6 个公开 bigWig 与 28 个样本的对应关系；deepTools 的归一化方法；公开 counts 的 peak 计数是否做过降采样。均登记 D7 与 gap 表。

---

## ④ 逐步可运行代码

> **实测数据格式**（2026-10-18 核对 GSE247685）：
> `GSE247685_consensus_peaks_kuramochi_all.tsv.gz` 不是纯 BED，而是 **chr/start/end + 20 个样本计数列**：
> `C_1,T5_1,T10_1,T20_1,T40_1,T80_1,T160_1,T320_1,C_2,T5_2,T10_2,T20_2,T40_2,T80_2,T160_2,T320_2,P10_1,P10_2,P320_1,P320_2`。
> 因此 **Fig. 3c–e 可以直接用公开数据严格复现**；缺的是 T320/T320W 的加药/撤药 8 个样本（Fig. 3f），需从 SRA 重跑。

### 4.1 主路线：公开 consensus counts → DESeq2

```r
# code/06_01_load_counts.R
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(DESeq2)
})
set.seed(20261018)

cnt <- fread("raw/geo/GSE247685/GSE247685_consensus_peaks_kuramochi_all.tsv.gz",
             data.table = FALSE, check.names = FALSE)
colnames(cnt) <- gsub("'", "", colnames(cnt))
stopifnot(all(c("chr","start","end") %in% colnames(cnt)))

peaks <- cnt[, c("chr","start","end")]
mat <- as.matrix(cnt[, setdiff(colnames(cnt), c("chr","start","end"))])
rownames(mat) <- paste0(peaks$chr, ":", peaks$start, "-", peaks$end)
storage.mode(mat) <- "integer"

meta <- data.frame(
  sample = colnames(mat),
  rep    = sub(".*_", "", colnames(mat)),
  condition = sub("_[12]$", "", colnames(mat)),
  stringsAsFactors = FALSE
)
print(table(meta$condition))

# Fig. 3c–e：C vs T5–T320（16 样本）；P10/P320 单独用于 persister 对照
keep <- meta$condition %in% c("C","T5","T10","T20","T40","T80","T160","T320")
dds <- DESeqDataSetFromMatrix(mat[, keep], meta[keep, ], design = ~ condition)
dds$condition <- relevel(factor(dds$condition), ref = "C")
dds <- DESeq(dds)

# 每个适应线 vs C
res_all <- lapply(setdiff(levels(dds$condition), "C"), function(cond) {
  r <- as.data.frame(results(dds, contrast = c("condition", cond, "C")))
  r$peak_id <- rownames(r); r$contrast <- paste0(cond, "_vs_C")
  cbind(peaks[match(r$peak_id, rownames(mat)), ], r)
}) %>% bind_rows()

res_all <- res_all %>%
  filter(!is.na(padj), padj < 0.001, abs(log2FoldChange) >= log2(3)) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "more", "less"))
write.csv(res_all, "results/06_diff_peaks.csv", row.names = FALSE)
cat("T320 less/more:",
    sum(res_all$contrast == "T320_vs_C" & res_all$direction == "less"), "/",
    sum(res_all$contrast == "T320_vs_C" & res_all$direction == "more"), "\n")
```

### 4.2 PCA（Fig. 3c）

```r
# code/06_02_pca.R
all_keep <- meta$condition %in% c("C","T5","T10","T20","T40","T80","T160","T320","P10","P320")
vsd <- vst(DESeqDataSetFromMatrix(mat[, all_keep], meta[all_keep, ], design = ~ condition))
top <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 2000)
pca <- prcomp(t(assay(vsd)[top, ]), scale. = TRUE)

pd <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                 condition = meta$condition[all_keep], sample = meta$sample[all_keep])
p <- ggplot(pd, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 2) +
  xlab(sprintf("PC1: %.1f%% variance", 100 * summary(pca)$importance[2,1])) +
  ylab(sprintf("PC2: %.1f%% variance", 100 * summary(pca)$importance[2,2])) +
  theme_classic()
ggsave("figures/Fig3c_atac_pca.pdf", p, width = 5, height = 4)
ggsave("figures/Fig3c_atac_pca.png", p, width = 5, height = 4, dpi = 300)
```

### 4.3 差异 peak 热图（Fig. 3d）

```r
# code/06_03_heatmap.R
library(ComplexHeatmap)
tp <- res_all %>% filter(contrast == "T320_vs_C") %>% pull(peak_id)
z <- assay(vsd)[tp, meta$condition[all_keep] %in% c("C","T5","T10","T20","T40","T80","T160","T320")]
z <- t(scale(t(z)))
ha <- HeatmapAnnotation(condition = meta$condition[all_keep][meta$condition[all_keep] %in% c("C","T5","T10","T20","T40","T80","T160","T320")])
pdf("figures/Fig3d_diff_peaks_heatmap.pdf", width = 8, height = 10)
Heatmap(z, top_annotation = ha, show_row_names = FALSE, name = "z-score")
dev.off()
```

### 4.4 motif 富集（Fig. 3e）

```bash
#!/usr/bin/env bash
# code/06_04_motif.sh
set -euo pipefail
OUT=results/06_motif
mkdir -p "$OUT"

# 从差异 peak 表生成 BED（T320 vs C）
Rscript -e '
r <- read.csv("results/06_diff_peaks.csv")
r <- r[r$contrast == "T320_vs_C", ]
for (d in c("more","less")) {
  x <- r[r$direction == d, c("chr","start","end")]
  write.table(x, paste0("results/06_motif/peaks_", d, ".bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
# Fig. 3f 用的合并 BED
write.table(r[, c("chr","start","end")], "results/06_motif/peaks_T320_vs_C.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
'

for d in more less; do
  findMotifsGenome.pl "results/06_motif/peaks_${d}.bed" hg38     "$OUT/${d}" -size 200 -mask -p 8
done

# 解析 HOMER knownResults.txt → 统一格式
Rscript -e '
library(dplyr); library(tidyr)
parse <- function(dir, direction) {
  f <- file.path(dir, "knownResults.txt")
  x <- read.delim(f, check.names = FALSE)
  data.frame(motif = x[["Motif Name"]], direction = direction,
             n_peaks = x[["# of Target Sequences with Motif"]],
             p = x[["P-value"]], p_adj = x[["q-value (Benjamini)"]])
}
bind_rows(parse("results/06_motif/more","more"),
          parse("results/06_motif/less","less")) %>%
  write.csv("results/06_motif_enrichment.csv", row.names = FALSE)
'
```

> HOMER 的列名随版本略有差异；若 `knownResults.txt` 列名不同，运行时打印 `colnames()` 后映射，不要静默改列。

### 4.5 CNA 控制（Extended Data Fig. 8c,d）

```r
# code/06_05_cna_control.R
library(GenomicRanges); library(dplyr)
res <- read.csv("results/06_diff_peaks.csv")
cna <- read.csv("results/05_cna_segments.csv")

gr <- GRanges(as.character(res$chr), IRanges(res$start, res$end))
cna_gr <- GRanges(as.character(cna$chr), IRanges(cna$start, cna$end))
hit <- countOverlaps(gr, cna_gr) > 0
res$overlaps_cna <- hit

# 分别生成 "no-CNA" 的 more/less BED，重跑 4.4 的 HOMER
for (d in c("more","less")) {
  x <- res[res$contrast == "T320_vs_C" & res$direction == d & !res$overlaps_cna,
           c("chr","start","end")]
  write.table(x, paste0("results/06_motif/peaks_", d, "_noCNA.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
cat("no-CNA fraction (less/more): ",
    mean(!res$overlaps_cna[res$direction == "less"]),
    mean(!res$overlaps_cna[res$direction == "more"]), "\n")
```

### 4.6 T320W 稳定性（Fig. 3f：需 SRA 重跑 8 个样本）

> 公开 consensus counts **不含** `T320_nodrug/T320_drug/T320W_nodrug/T320W_drug`。若要复现 Fig. 3f，必须从 **SRA/BioProject PRJNA1040167** 下载这 8 个样本的 fastq，按 4.7 的流程比对、call peak、计数。

```bash
#!/usr/bin/env bash
# code/06_06_t320w.sh
set -euo pipefail
OUT=processed/06_t320w
mkdir -p "$OUT"/{bam,bigwig}
# 前提：8 个样本已按 06_07 流程生成 dedup.bam
# 样本映射：T320 = T320_nodrug（撤药 10 天）；T320W = T320W_nodrug（撤药 120 天）
for s in C T320_nodrug T320_drug T320W_nodrug T320W_drug; do
  bamCoverage -b "$OUT/bam/${s}.dedup.bam" -o "$OUT/bigwig/${s}.bw" \
    --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p 12
done

# 用 Fig. 3d 的 T320 vs C 差异 peak（±1 kb）做热图
computeMatrix reference-point \
  -S "$OUT"/bigwig/{C,T320_nodrug,T320_drug,T320W_nodrug,T320W_drug}.bw \
  -R results/06_motif/peaks_T320_vs_C.bed \
  --referencePoint center -b 1000 -a 1000 --binSize 1 \
  -o "$OUT/t320w_matrix.gz" -p 12

plotHeatmap -m "$OUT/t320w_matrix.gz" -o figures/Fig3f_T320W_heatmap.pdf \
  --colorMap RdBu_r --whatToShow 'heatmap and colorbar'
plotHeatmap -m "$OUT/t320w_matrix.gz" -o figures/Fig3f_T320W_heatmap.png \
  --colorMap RdBu_r --whatToShow 'heatmap and colorbar' --dpi 300
```

```r
# Wilcoxon：less/more peak 在 T320 vs T320W 的信号差
library(rtracklayer); library(GenomicRanges)
sig <- function(bw, gr) {
  x <- import(bw, which = gr, as = "NumericList")
  sapply(x, mean, na.rm = TRUE)
}
# 对 more peaks 与 less peaks 分别比较；报告 Cohen's d
```

### 4.7 仅需重跑时的 SRA 流程（Fig. 3f / 需要原始比对）

```bash
#!/usr/bin/env bash
# code/06_07_atac_pipeline.sh
set -euo pipefail
REF=raw/reference
OUT=processed/06_atac
mkdir -p "$OUT"/{trim,bam,peaks,bigwig,logs}

for s in $(ls raw/fastq/*_R1.fastq.gz | xargs -n1 basename | sed 's/_R1.fastq.gz//'); do
  trim_galore --paired -q 20 --length 30 --gzip \
    -o "$OUT/trim" "raw/fastq/${s}_R1.fastq.gz" "raw/fastq/${s}_R2.fastq.gz"

  bowtie2 --very-sensitive -k 10 -X 2000 -p 12 -x "$REF/GRCh38" \
    -1 "$OUT/trim/${s}_R1_val_1.fq.gz" -2 "$OUT/trim/${s}_R2_val_2.fq.gz" \
    2> "$OUT/logs/${s}.bowtie2.log" \
  | samtools view -@ 8 -bS - \
  | samtools sort -@ 8 -o "$OUT/bam/${s}.sorted.bam" -

  samtools index "$OUT/bam/${s}.sorted.bam"
  samtools idxstats "$OUT/bam/${s}.sorted.bam" | cut -f1 | grep -v '^chrM$' > "$OUT/logs/${s}.chroms"
  samtools view -@ 8 -b "$OUT/bam/${s}.sorted.bam" $(cat "$OUT/logs/${s}.chroms" | tr '\n' ' ') \
    > "$OUT/bam/${s}.noMT.bam"

  picard MarkDuplicates I="$OUT/bam/${s}.noMT.bam" O="$OUT/bam/${s}.dedup.bam" \
    M="$OUT/logs/${s}.dup_metrics.txt" REMOVE_DUPLICATES=true
  samtools index "$OUT/bam/${s}.dedup.bam"
done

# 记录深度并降采样到最低深度（samtools view -s）
for bam in "$OUT"/bam/*.dedup.bam; do
  s=$(basename "$bam" .dedup.bam)
  samtools view -c -F 1804 "$bam" | awk -v s="$s" '{print s"\t"$1}'
done > "$OUT/logs/library_depth.tsv"
MIN=$(cut -f2 "$OUT/logs/library_depth.tsv" | sort -n | head -1)
while IFS=$'\t' read -r s n; do
  frac=$(awk -v n="$n" -v m="$MIN" 'BEGIN{printf "%.6f", m/n}')
  samtools view -b -s "$frac" "$OUT/bam/${s}.dedup.bam" > "$OUT/bam/${s}.down.bam"
  samtools index "$OUT/bam/${s}.down.bam"
done < "$OUT/logs/library_depth.tsv"

# MACS2 + IDR + consensus + deepTools 见 01 号文档的参考流程
```

---

## ⑤ 参数解释与调参空间

| 参数 | 原文 | 本项目 | 改变后果 |
|---|---|---|---|
| 细胞数/条件 | 50,000 | 50,000 | 过少导致 peak 稀疏 |
| TrimGalore | q20, length 30 | 同 | 过严丢失短读 |
| Bowtie2 | very-sensitive, -k 10, -X 2000 | 同 | 默认更慢但更敏感 |
| 降采样 | 最低深度 | 同 | 不降采样会引入深度偏差 |
| MACS2 q | 0.05 | 0.05 | 影响 peak 数量 |
| IDR | <0.05 | 同 | 过松会保留假 peak |
| DESeq2 阈值 | adj P<0.001, FC≥3 | 同 | 过松增加假阳性 |
| PCA peak 数 | 2,000 | 2,000 | 影响方差解释比例 |
| motif 背景 | HOMER 默认 | 同 | 不同背景会改变富集 |
| T320W 撤药 | 120 天 | 120 天 | 时间越长，可逆性越明显 |
| 置换 | 1,000 | 1,000 | CNA 重叠检验 |
| 归一化 | RPGC | 同 | 影响跨样本比较 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 原文报告 | 判据 |
|---|---|---|
| 样本数 | 公开 counts 20 个（10 条件 × 2 重复）；全实验 28 个 | 公开计数表列名核对；缺 8 个 T320/T320W 加药样本 |
| PCA PC1/PC2 | 80% / 8% | PC1 沿 C→T320 单调；方差 ±10 个百分点 |
| 差异 peak（T320 vs C） | less 6,771；more 7,332 | 同数量级（±30%） |
| motif 下降 | SOX17、PAX8、WT1 | 方向一致 |
| motif 上升 | AP1、NRF2、ATF4 | 方向一致 |
| CNA 控制 | 70%（less）/80–90%（more）不重叠 CNA | 重算后方向不变 |
| T320W | 10/40 μM 耐受、320 μM 敏感 | 趋势一致 |
| 再暴露 | 恢复开放 | T320W_drug > T320W_nodrug |
| 增殖 | T320/T320W 不慢于 C | 无适应性代价 |
| 重复性 | 两次运行 | peak 集合 Jaccard >0.8 |

**重要限制**：公开 counts 覆盖 20/28 样本。Fig. 3c–e 可严格复现；**Fig. 3f（T320W 加药/撤药）无法用公开数据复现**，必须走 SRA。gap 表必须写明这一点，不要把"没做"写成"未复现"。

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 找不到 T320W 的计数 | 公开 counts 只有 20 样本，缺 8 个加药/撤药样本 | 走 SRA（D7）；不要用 bigWig 冒充计数 |
| 把 6 个 bigWig 当 28 样本 | 系列级轨道，无法按样本拆分 | 只用 consensus counts 做 DE |
| 线粒体读段占比高 | 核提取不纯 | 加 digitonin；去 MT reads |
| peak 数与原文差很多 | MACS2 参数/深度不同 | 固定 q、nomodel、extsize；报告差异 |
| IDR 报错 | narrowPeak 排序不一致 | 用 `--rank p.value`；确保同条件重复配对 |
| DESeq2 报非整数 | deepTools 输出浮点 | `round()`；检查是否用了 bigWig 信号 |
| PCA 不分条件 | 批次/深度主导 | 检查降采样与归一化 |
| motif 富集全是垃圾 | 背景选择不当 | 用 `-bg` 指定匹配背景；`-mask` |
| CNA 控制后结果反转 | CNA 是主要驱动 | 如实报告；不要隐藏 |
| T320W 信号不恢复 | 撤药/再暴露时间不够 | 核对 120 天与 320 μM 处理 |
| 内存不足 | 1 bp bin 的 bigWig | 增大 binSize 或在目标区域计算 |

---

## ⑧ 与文献图表的对应关系

| 原文图 | 内容 | 本项目产物 | 说明 |
|---|---|---|---|
| Fig. 3b | ATAC 实验设计（分选） | `figures/Fig3b_atac_scheme.pdf/png` | 方案图 |
| Fig. 3c | PCA（PC1 80%、PC2 8%） | `figures/Fig3c_atac_pca.pdf/png` | top 2,000 peak |
| Fig. 3d | 差异 peak 热图（T5–T320 vs C） | `figures/Fig3d_diff_peaks_heatmap.pdf/png` | 路线 A |
| Fig. 3e | motif 富集条形图 | `figures/Fig3e_motif_bars.pdf/png` | SOX17/PAX8/WT1↓；AP1/NRF2/ATF4↑ |
| Fig. 3f | T320 vs T320W ± 药热图 | `figures/Fig3f_T320W_heatmap.pdf/png` | 1 bp、RPGC |
| Extended Data Fig. 8a | 差异 peak vs 总 peak 相关 | `results/06_diff_vs_total_peaks.csv` | 排除检测偏差 |
| Extended Data Fig. 8b | motif 热图 | `figures/FigED8b_motif_heatmap.pdf/png` | |
| Extended Data Fig. 8c | peak–CNA 重叠 | `results/06_peak_cna_overlap.csv` | 1,000 置换 |
| Extended Data Fig. 8d | 去 CNA 后 motif | `figures/FigED8d_motif_noCNA.pdf/png` | 关键对照 |
| Extended Data Fig. 8e | T320W 增殖/药敏 | `figures/FigED8e_T320W_growth.pdf/png` | 09 号模块 |
| Extended Data Fig. 8f | less/more peak 信号箱线图 | `figures/FigED8f_peak_signal.pdf/png` | Wilcoxon、Cohen's d |

> **交付判定**：用公开 counts 需达到 PCA 顺序、motif 方向、差异 peak 数量级三项（Fig. 3c–e）；T320W（Fig. 3f）只在完成 SRA 重跑后判定，否则在 10 号文档 gap 表标注"数据未公开，未执行"。
