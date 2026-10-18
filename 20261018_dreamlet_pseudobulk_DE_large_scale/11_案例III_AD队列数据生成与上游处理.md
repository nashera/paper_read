# 11 案例 III：AD 队列数据生成与上游处理

> **目标文献**：Hoffman GE, et al. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* Nature Communications (2026) 17:9597. DOI: 10.1038/s41467-026-75680-8
> **模块对应原文**：Fig. 5A–B；Supplementary Fig. 6；Methods「Single-nucleus RNA-seq data generation」全部小节（Study cohort、Isolation and FANS、snRNA-seq and library preparation、Single-nucleus data processing and quality control）
> **前置依赖**：`01_环境搭建与数据获取.md`（`raw/psychad_ad/`）；若只做下游，可跳过本模块直接用作者 h5ad
> **复现难度 / 预计工时 / 计算资源**：★★★★★（湿实验）/ ★★★★☆（干实验）/ 湿实验数周–数月；干实验 2–5 天 / 集群数百核时

---

## ① 科学问题与逻辑位置

- **上游**：尸检 DLPFC 脑组织（Mount Sinai NIH NeuroBioBank，>2000 例脑）。
- **下游**：`processed/sce_psychad_raw.rds`（细胞 × 基因计数 + 22 类注释 + 元数据）→ 12 号。
- **本模块回答**：**"1.4 M 个单核转录组是怎么从冰冻脑组织一步步变成可用于跨受试者 DE 的注释对象的？每一步的偏差在哪里？"**

> ⚠️ **导师批注**：这是全套文档里**唯一无法在普通实验室完整重跑**的模块（需要尸检脑组织、FANS 设备、10x、NovaSeq、集群）。本项目的默认路径是**使用作者公开的 h5ad**（D1），本模块的价值是：① 让你理解数据是怎么来的；② 在需要"从 FASTQ 重跑"时有一份可执行的 SOP；③ 让你在审稿/质疑时能指出上游每一步的自由度。

---

## ② 背景原理

### 2.1 为什么用单核（snRNA-seq）而不是单细胞

| 维度 | scRNA-seq | snRNA-seq |
|---|---|---|
| 样本类型 | 新鲜/冻存细胞 | **冰冻组织**（尸检脑只能冻存） |
| 分离难度 | 需要活细胞、酶解 | 只需裂解细胞膜，核完整即可 |
| 偏倚 | 酶解偏好、细胞大小 | 核 RNA 偏核内（内含子、lncRNA） |
| 适用 | 血液、类器官 | **脑、脂肪、肌肉**等难解离组织 |

冰冻尸检脑无法用酶解获得高质量单细胞，所以原文用**蔗糖梯度 + Triton X-100 裂解**提取核。

### 2.2 多重混合（multiplexing）与 hashing

**问题**：586 个样本逐个建库，成本高、批次多。

**方案**：6 个供体一组，每个样本用**不同的 TotalSeq-A 抗核抗体**（带独特 barcode 的 oligo）标记，混合后一起建库。

```text
样本1 + hash-A ─┐
样本2 + hash-B ─┤
样本3 + hash-C ─┼─→ 混合 60,000 核 ─→ 10x 3' v3.1 ─→ 一个文库
样本4 + hash-D ─┤
样本5 + hash-E ─┤
样本6 + hash-F ─┘
```

**为什么还要基因型去多重？** Hashing 抗体可能脱落、交叉反应；原文同时用 **cellSNP + vireo** 从表达数据中的 SNP 推断每个细胞属于哪个 donor，再用 **QTLtools-mbv** 与 SNP array 基因型比对，纠正样本错标。这是"双保险"。

### 2.3 关键批次结构

| 层级 | 变量 | 数量级 | 处理 |
|---|---|---|---|
| 供体 | `SubID` | 299 | 随机效应 |
| 10x pool（6 供体混合） | `poolID` | 以实际列为准（D18） | 随机效应（高维批次） |
| 建库通道 | `Channel` | 586 | **样本单位**（pseudobulk 的 `sample_id`） |
| 测序轮次 | `round_num` | 数个 | 记录，必要时随机效应 |

> ⚠️ **导师批注**：`Channel`（一次 10x 建库）≠ `SubID`（人）。同一供体可能有多个 Channel（技术重复）。这正是原文 Fig. 5G "technical replicates" 的来源：同一 SubID 的不同 Channel 之间的一致性，就是"技术可重复性"。

### 2.4  Ambient RNA、双胞与空滴

- **Ambient RNA**：裂解时释放的游离 mRNA 被其他核捕获 → 低表达基因被"污染"。原文检查 rRNA、sRNA、假基因、MALAT1 的比例。
- **Doublets**：两个核被同一个 barcode 捕获 → 表达谱是混合。Scrublet 用模拟双胞的 kNN 距离检测。
- **空滴**：只有 ambient RNA、没有核 → UMI 很低，QC 过滤。

---

## ③ 原文逐句拆解

### 3.1 队列纳入标准

| 原文句 | 中文转述 | 参数 |
|---|---|---|
| "AD samples needed to be classified by (1) CERAD protocol as 'AD possible', 'AD probable' or 'AD definite', (2) Braak AD staging protocol within the stages 3–6, and (3) clinical dementia rating within the rate 0.5-5." | AD 三条件 | CERAD / Braak 3–6 / CDR 0.5–5 |
| "our subset of AD donors cannot be diagnosed with Parkinson's disease or Diffuse Lewy body disease." | 排除 PD/DLB | — |
| "control samples needed to be classified by (1) CERAD protocol as 'no AD' or 'possible AD' and (2) Braak AD staging protocol within the stage 0-2." | 对照两条件 | CERAD no/possible；Braak 0–2 |
| "control samples cannot be diagnosed with any other neurodegenerative, neurological or neuropathological diagnosis." | 对照排除其他神经疾病 | — |
| "donors over 60 years of age" | 年龄 >60 | — |

### 3.2 湿实验

| 原文句 | 中文转述 | 关键参数 |
|---|---|---|
| "25 mg of frozen postmortem human brain tissue … homogenized in cold lysis buffer (0.32 M Sucrose, 5 mM CaCl2, 3 mM Magnesium acetate, 0.1 mM EDTA, 10 mM Tris-HCl, pH8, 1 mM DTT, 0.1% Triton X-100) and filtered through a 40 µm cell strainer." | 裂解液配方 | 见左 |
| "underlaid with sucrose solution (1.8 M Sucrose, 3 mM Magnesium acetate, 1 mM DTT, 10 mM Tris-HCl, pH8) and centrifuged at 107,000 g for 1 hour at 4 °C." | 蔗糖 cushion 超速离心 | 107,000 g / 1 h |
| "2 M from each sample were pelleted at 500 g for 5 min at 4 °C and re-suspended in 100 µl staining buffer (2% BSA, 0.02% Tween-20 in PBS)." | 每样本 2 M 核 | — |
| "Each sample incubated with 1 µg of a distinct TotalSeq-A nuclear hashing antibody (Biolegend) for 30 min at 4 °C." | hashing | 1 µg / 30 min |
| "7-AAD positive nuclei were collected … FACSAria" | FANS 分选 | 7-AAD⁺ |
| "Two aliquots of 60,000 pooled nuclei (i.e., 10,000 each) were processed in parallel using 3' v3.1 reagents." | 两个 60k 核的通道 | 10,000/样本 |
| "At the cDNA amplification step (step 2.2), reactions were supplemented with a hash-tag oligo (HTO) cDNA 'additive' primer (GTGACTGGAGTTCA-GACGTGTGCTCTTCCGAT*C*T; *Phosphorothioate bond)." | HTO 引物 | 序列见左 |
| "supernatants from the 0.6x SPRI selection step were retained for HTO library generation." | HTO 文库 | 0.6× SPRI |
| "All libraries were sequenced at NYGC using the Novaseq platform." | 测序 | NovaSeq |

### 3.3 干实验

| 原文句 | 中文转述 | 参数 |
|---|---|---|
| "Sequencing reads … aligned to the hg38 reference genome using STARsolo." | STARsolo | hg38 |
| "cellSNP was used to pile up the expressed alleles from polymorphic sites overlapping snRNA-seq reads." | cellSNP | 未给阈值 |
| "vireo utilized those pile-ups to split cells into clusters corresponding to six distinct donors per pool." | vireo，每 pool 6 donor | `-N 6` |
| "genotype concordance analysis that compared the clusters of cells against reference SNP-array data using QTLtools-mbv." | QTLtools-mbv | — |
| "Poor-quality cells were detected by thresholding based on UMI counts, gene counts, and mitochondrial contents." | 三级 QC | **未给具体阈值** |
| "removing doublets using the Scrublet method." | Scrublet | 未给阈值 |
| "We removed features (genes) that are not robustly expressed by at least 0.05% of the cells/nuclei." | 基因过滤 | 0.05% |
| "We remove donors with fewer than 50 cells." | 供体过滤 | ≥50 |
| "Harmony … to regress out the effect of sequencing pools before performing clustering." | Harmony | pool |
| "kNN graph … Leiden clustering … UMAP" | 聚类 | 未给分辨率 |
| "annotated based on a combination of expert curation and machine-learning-based algorithms to query known gene marker signatures previously curated by Human Cell Atlas." | 注释 | — |
| "For each donor, the labeled sex was checked to be consistent with expression on XIST and UTY." | 性别核查 | XIST/UTY |

---

## ④ 逐步可运行代码

### 4.1 STARsolo 比对与定量

```bash
# code/11_01_starsolo.sh
# 10x 3' v3.1：CB 16 bp，UMI 12 bp；R1 是 barcode+UMI，R2 是 cDNA
STAR --runThreadN 16 \
  --genomeDir ref/STAR_hg38 \
  --readFilesIn fastq/${SAMPLE}_R2.fastq.gz fastq/${SAMPLE}_R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist ref/3M-february-2018.txt \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ \
  --soloMultiMappers EM \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix out/${SAMPLE}_
```

> `--soloFeatures Gene GeneFull` 同时输出"外显子计数"与"全基因（含内含子）计数"。snRNA-seq 常用 `GeneFull`，但**作者代码读的是 `Solo.out/Gene/raw`**，故本项目主路径用 **Gene**（D19）；`GeneFull` 只作敏感性分析，并在 QC 表记录两者基因数差异。

### 4.2 基因型去多重

```bash
# code/11_02_demux.sh
# (1) cellSNP：在表达数据中找 SNP
cellsnp-lite -s out/${SAMPLE}_Aligned.sortedByCoord.out.bam \
  -b out/${SAMPLE}_Solo.out/Gene/raw/barcodes.tsv \
  -O cellSNP/${SAMPLE} -p 16 --minMAF 0.1 --minCOUNT 20

# (2) vireo：每 pool 6 个供体
vireo -c cellSNP/${SAMPLE} -N 6 -o vireo/${SAMPLE} --randSeed 20261018

# (3) QTLtools-mbv：与 SNP array 比对，纠正错标
QTLtools mbv --bam out/${SAMPLE}_Aligned.sortedByCoord.out.bam \
  --vcf genotype/${SAMPLE}.vcf.gz --out mbv/${SAMPLE}
```

> ⚠️ **导师批注**：`-N 6` 是**假设**每 pool 恰好 6 个供体（原文设计如此）。若某 pool 实际只有 5 个或混入 7 个，vireo 会强行分 6 类，产生"幽灵供体"。必须用 QTLtools-mbv 的一致性结果反查。

### 4.3 R 端 QC（三级）

```r
# code/11_03_qc.R
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(scater); library(scuttle)
  library(scran); library(scDblFinder); library(BiocParallel)
})

sce <- readRDS("processed/sce_psychad_raw.rds")   # 已合并所有 Channel

# ★ 基因 ID 统一：h5ad 是 Ensembl（D2），先按 id_map / rowData 转 SYMBOL
if (grepl("^ENSG", rownames(sce)[1])) {
  sym <- rowData(sce)$gene_name
  if (is.null(sym)) {
    map <- read.delim("processed/id_map_psychad.tsv")
    sym <- map$symbol[match(rownames(sce), map$original_id)]
  }
  sym[is.na(sym) | sym == ""] <- rownames(sce)[is.na(sym) | sym == ""]
  rownames(sce) <- scater::uniquifyFeatureNames(rownames(sce), sym)
}
stopifnot("MALAT1" %in% rownames(sce), any(grepl("^MT-", rownames(sce))),
          "XIST" %in% rownames(sce), "UTY" %in% rownames(sce))

# ---- 细胞级 ----
# 显式定义 MT 子集（此时行名已是 SYMBOL）
is_mito <- grepl("^MT-", rownames(sce))
qc <- addPerCellQCMetrics(sce, subsets = list(MT = is_mito))

# ⚠️ 原文未给具体阈值（D16）；以下为本项目默认，必须做敏感性分析
umi_min   <- 500
gene_min  <- 250
mito_max  <- 0.05
sce$pass_cell <- qc$sum > umi_min &
                 qc$detected > gene_min &
                 qc$subsets_MT_percent < 100 * mito_max

# ambient RNA 代理指标：MALAT1 占比（核内高表达，过高提示 ambient）
sce$frac_MALAT1 <- as.numeric(counts(sce)["MALAT1", ]) / pmax(qc$sum, 1)

# ---- 双胞 ----
# 原文用 Scrublet；R 侧用 scDblFinder 作为替代（不同算法，需在方法学声明）
dbl <- scDblFinder::scDblFinder(sce, samples = "Channel", BPPARAM = BiocParallel::SnowParam(8))
sce$doublet <- dbl$scDblFinder.class
sce <- sce[, sce$doublet == "singlet"]

# ---- 基因级 ----
keep_gene <- rowSums(counts(sce) > 0) >= 0.0005 * ncol(sce)   # 0.05%
sce <- sce[keep_gene, ]

# ---- 供体级 ----
tab <- table(sce$SubID)
keep_donor <- names(tab)[tab >= 50]
sce <- sce[, sce$SubID %in% keep_donor]

# ★ 保存 QC 后对象（12/13 号的输入，契约表登记）
saveRDS(sce, "processed/sce_psychad_qc.rds")
writeLines(capture.output(sessionInfo()), "logs/11_qc_sessionInfo.txt")
```

> **必须登记**：UMI/基因/线粒体阈值是**复现者决策**（D16）：原文对这些只写 "thresholding"，没给数字；但基因 0.05% 与供体 ≥50 细胞是原文明确写的，按原文执行。敏感性分析：把 `umi_min` 改为 300/1000，看 DEG 数与 PTPRG 效应是否稳定。

### 4.4 Harmony + Leiden + UMAP

```r
# code/11_04_cluster.R
suppressPackageStartupMessages({
  library(scran); library(scater); library(harmony); library(igraph); library(leiden)
})

sce <- readRDS("processed/sce_psychad_qc.rds")

# HVG（按均值–方差趋势；与原文一致）
dec <- modelGeneVarByPoisson(sce)
hvg <- getTopHVGs(dec, n = 2000)

# PCA → Harmony（按 pool）
sce <- runPCA(sce, subset_row = hvg, ncomponents = 30)
sce <- harmony::RunHarmony(sce, group.by.vars = "poolID",
                           reduction = "PCA", reduction.save = "HARMONY",
                           seed = 20261018)

# kNN 图 + Leiden（bluster 封装，避免直接依赖 igraph 的版本差异）
cl <- bluster::clusterRows(
  reducedDim(sce, "HARMONY"),
  bluster::NNGraphParam(k = 20, cluster.fun = "leiden", cluster.args = list(resolution = 1))
)
sce$cluster <- factor(cl$clusters)
```

> `bluster` 的 `cluster.fun = "leiden"` 需要 `leiden` 包（`install.packages("leiden")`）。若不可用，退化为 `"louvain"`（`igraph` 内置），并在 `DECISIONS.md` 记录。**所有聚类结果都必须在 `logs/` 记录分辨率与 k**。

### 4.5 注释与性别核查

```r
# code/11_05_annotate.R
suppressPackageStartupMessages({ library(scran); library(scater); library(SingleCellExperiment) })
sce <- readRDS("processed/sce_psychad_qc.rds")

# 主路径：public release 已提供 subclass（22 类），直接沿用（D12）
if ("subclass" %in% colnames(colData(sce))) {
  sce$subclass <- factor(sce$subclass)
  message("using provided subclass labels: ", length(levels(sce$subclass)))
} else {
  # 备选：从 FASTQ 重跑时，用 marker 打分做**粗注释**（10 大类），再人工细分 22 亚型
  sce <- logNormCounts(sce)
  markers <- list(
    Micro_PVM = c("CSF1R","CX3CR1","P2RY12","TMEM119"),
    Astro     = c("GFAP","AQP4","SLC1A2"),
    Oligo     = c("PLP1","MBP","MAG"),
    OPC       = c("PDGFRA","VCAN","OLIG2"),
    Endo      = c("FLT1","CLDN5","PECAM1"),
    VLMC      = c("PDGFRB","SLC14A1"),
    PC        = c("PDGFRB","RGS5"),
    CD8_T     = c("CD3D","CD8A","CD8B"),
    EN        = c("SLC17A7","RBFOX3"),
    IN        = c("GAD1","GAD2")
  )
  sm <- scoreMarkers(sce, groups = sce$cluster)
  ct_score <- sapply(names(markers), function(ct) {
    g <- intersect(markers[[ct]], rownames(sce))
    vapply(names(sm), function(cl) mean(sm[[cl]][g, "mean.AUC"], na.rm = TRUE), numeric(1))
  })
  coarse <- setNames(colnames(ct_score)[max.col(ct_score, ties.method = "first")],
                     rownames(ct_score))
  sce$coarse <- coarse[as.character(sce$cluster)]
  message("coarse annotation done; 22 subtypes require expert curation")
}

saveRDS(sce, "processed/sce_psychad_annot.rds")
writeLines(capture.output(sessionInfo()), "logs/11_annot_sessionInfo.txt")
```

```r
# code/11_06_sex_check.R
# 按 SubID 汇总（pseudobulk 求和），而不是对细胞取等权均值
sce <- readRDS("processed/sce_psychad_annot.rds")
xist <- as.numeric(counts(sce)["XIST", ])
uty  <- as.numeric(counts(sce)["UTY", ])
lib  <- colSums(counts(sce))
df <- data.frame(SubID = sce$SubID, sex_label = sce$Sex,
                 XIST = xist, UTY = uty, lib = lib)
agg <- aggregate(cbind(XIST, UTY, lib) ~ SubID + sex_label, data = df, FUN = sum)
agg$XIST_rate <- agg$XIST / pmax(agg$lib, 1)
agg$UTY_rate  <- agg$UTY  / pmax(agg$lib, 1)
agg$call <- ifelse(agg$XIST_rate > agg$UTY_rate, "F", "M")
print(table(agg$sex_label, agg$call))
write.csv(agg, "results/11_sex_check.csv", row.names = FALSE)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| 组织量 | 25 mg | 10–50 mg | 核产量与 ambient RNA |
| 蔗糖浓度 | 1.8 M | 1.5–2.0 M | 核纯度 |
| 离心 | 107,000 g / 1 h | — | 标准 |
| 每样本核数 | 2 M | — | hashing 效率 |
| 混合规模 | 6 样本 / 60k 核 | 4–12 | pool 数与批次效应 |
| STARsolo `soloFeatures` | Gene | GeneFull | 内含子计数比例 |
| vireo `-N` | 6 | 实际供体数 | 幽灵供体 |
| UMI 阈值 | **500（我方）** | 300–1000 | 细胞数/DEG 数（D16） |
| 基因阈值 | **250（我方）** | 200–500 | 同上 |
| 线粒体 | **5%（我方）** | 5–10% | 核质量 |
| 基因表达比例 | 0.05% | 0.01–0.1% | 基因数 |
| 供体细胞数 | 50 | 20–100 | 供体数 |
| Harmony 维度 | 30 | 20–50 | 聚类粒度 |
| Leiden resolution | 1.0 | 0.4–2.0 | 簇数（决定 22 类是否成立） |

> ⚠️ **导师批注**：原文的 QC 阈值、Leiden 分辨率、注释方法都只给了定性描述。**22 类**是原文的结果，不是参数。复现时若得到 20 或 25 类，不要强行调到 22——应先比对 marker 与已知 HCA 图谱，再决定是"合并/拆分"还是"参数不同"。这就是 D16 的敏感性分析要解决的问题。

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| 样本数 | 586 | ±5% | ❌ |
| 供体数 | 299 | ±5% | ❌ |
| 细胞核数 | 1.4 M | ±5% | ❌ |
| AD / Control | 150 / 149 | ±5% | ❌ |
| 细胞类型数 | 22（8 EN + 6 IN + 8 其他） | 必须可映射 | ❌ |
| 年龄 | >60 | 必须 | ❌ |
| 性别核查一致率 | >95% | — | ⚠️ |
| 双胞率 | 5–15%（视组织） | — | ⚠️ |
| MALAT1 比例 | 高（核内） | — | ✅ 预期 |
| 每供体细胞数中位 | 数千 | — | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| vireo 分出 7 类 | `-N` 设错或混入额外样本 | 用 `-N` 实际值；QTLtools-mbv 校验 |
| 某 Channel 全是 doublet | 核浓度过高 | 稀释后重跑；或 Scrublet 过滤 |
| 线粒体比例整体很高 | 核膜破裂 | 优化裂解时间/Triton 浓度 |
| 基因数很少 | 用了 Gene 而非 GeneFull | 换 `soloFeatures` |
| 聚类出现"按 pool 分"的簇 | 未做 Harmony 或 Harmony 未收敛 | 检查 `group.by.vars`；增加迭代 |
| 注释混乱 | marker 面板不匹配组织 | 用 HCA 脑图谱 marker；专家复核 |
| sex check 不一致 | 样本错标 | 回到 QTLtools-mbv 的基因型一致性 |
| 22 类对不上 | 分辨率/注释不同 | 先映射到 22 类，再比较；不要硬调 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 5A（实验设计示意） | `figures/Fig5A_design.*` | 6 供体/pool、两个技术重复 |
| Fig. 5B（UMAP + 22 类） | `figures/Fig5B_umap_22ct.*` | `plotProjection()` |
| Fig. S6（人口学、QC） | `figures/FigS6_demographics.*` | 年龄/性别/细胞数 |
| Methods 湿实验 | `logs/11_wetlab_protocol.md` | 本文件 §3.2 + §4.1 |
| Methods 干实验 | `logs/11_drylab_commands.sh` | §4.1–4.6 |

> **交付物**：`processed/sce_psychad_qc.rds`、`processed/sce_psychad_annot.rds`、`figures/Fig5A/B`、`logs/11_*`。完成后进入 `12_案例IV_AD差异表达与下游分析.md`。
