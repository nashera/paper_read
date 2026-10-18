# 06 细胞互作网络 I：scRNA-seq × TCGA/GTEx 相关性网络

> **目标文献**：Zhang L, Li Z, Skrzypczynska KM, et al. *Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer.* Cell. 2020;181(2):442–459.e29. DOI: 10.1016/j.cell.2020.03.048
> **模块对应原文**：Fig. 3E；Fig. S5A–S5C；Fig. S8A；STAR Methods「Correlative cell-cell interactions inferred by combined scRNA-seq and TCGA datasets」（p.e11）
> **前置依赖**：`processed/human_smartseq2_annot.rds`（03）、`processed/human_smartseq2_norm.rds`（02）、`results/ro_e.tsv`（03）、`results/myeloid_markers.tsv`（04）、`processed/human_10x_norm.rds`（02）；`raw/tcga_coadread/`、`raw/gtex_colon/`（01）；`PROJECT_SPEC.md` §4 契约、§5 术语、D11
> **复现难度 / 预计工时 / 计算资源**：★★★★☆ / 5–7 人日 / 32 GB 内存、8 CPU；TCGA/GTEx 相关矩阵可并行

---

## ① 科学问题与逻辑位置

### 1.1 这一模块要回答什么

**在只有 scRNA-seq（细胞类型、但样本量小）和 bulk RNA-seq（样本量大、但细胞类型混在一起）的情况下，能否推断"哪些细胞类型在肿瘤中共同出现/相互依赖"，并找出网络的核心节点？**

原文的回答：可以，但只能得到**相关性网络（correlative network）**，不是物理互作网络。做法是：

1. 用 scRNA-seq 定义每个细胞亚群的 **signature 基因**；
2. 把这些 signature 投射到 TCGA COAD/READ（肿瘤）与 GTEx colon（正常黏膜）的 bulk 表达，估计每个样本中每个亚群的**相对丰度**；
3. 对每个基因与每个亚群丰度做 **Pearson 相关**，得到"基因–细胞丰度"矩阵；
4. 去掉 **self-expressed genes**（自己定义自己的基因），取每个亚群相关最高的 **top 13** 非自表达基因；
5. 把这些基因放回 scRNA-seq，看它们在哪一个亚群中高表达（**富集分数 z>1.96**），从而连边；
6. 肿瘤与正常分别建网，比较网络结构。

### 1.2 输入 / 输出

| 类型 | 内容 | 契约 |
|---|---|---|
| 输入 | Smart-seq2 `logtpm`（基因 × 细胞）、`human_smartseq2_annot.rds` 的 `cluster_id`/`tissue` | `PROJECT_SPEC.md` §4 |
| 输入 | TCGA COAD+READ、GTEx colon 的 bulk 表达（gene symbol × sample，log2 尺度） | 01 号下载并校验 |
| 输入 | 每个簇的 Ro/e 组织富集表 | `results/ro_e.tsv`（03 号；若 `PROJECT_SPEC.md` §4 的实际契约名不同，以契约为准，不在脚本里猜路径） |
| 输出 | `results/celltype_signatures.tsv` | 06 号唯一生产者 |
| 输出 | `results/abundance_tcga.tsv`、`results/abundance_gtex.tsv` | 06 号唯一生产者 |
| 输出 | `results/correlation_gene_cluster.tsv`（长表）、`results/correlation_network_edges.tsv` | 06 号唯一生产者 |
| 输出 | `figures/Fig3E_network.pdf/.png`、`figures/FigS5C_network_normal.pdf/.png` | 06 号唯一生产者 |

### 1.3 在证据链中的位置

```text
02 预处理 → 03 聚类/注释/Ro/e → 04 髓系亚群/signature
                                      ↓
06 相关性网络（本模块） ← TCGA/GTEx bulk
                                      ↓
07 配体–受体（分子层） → 08 小鼠模型验证 → 12 生存验证
```

本模块是**连接"人类图谱"与"bulk 临床队列"的桥梁**：没有它，后续的 C1QC/SPP1 生存分析（09/12）与 LR 分析（07）就缺少"哪些细胞类型值得看"的依据。

---

## ② 背景原理

### 2.1 为什么用 bulk 而不是直接做细胞间互作

- scRNA-seq 的样本量通常只有几十例，无法做"细胞类型丰度 × 临床/分子特征"的稳健相关；
- TCGA/GTEx 有数百例，但每个样本是混合细胞，无法直接看到细胞类型；
- **折中方案**：用 scRNA-seq 学到"每个细胞类型的基因指纹"，再用这个指纹在 bulk 中反推"相对丰度"。这不是严格的去卷积（deconvolution），而是**signature score**：假设 signature 基因的平均表达与细胞丰度单调相关。
- 因此，得到的网络是**共变网络**：A 丰度高时 B 的 marker 也高，可能是 A–B 互作，也可能是两者被同一因素驱动（肿瘤纯度、炎症状态、样本构成）。原文用 "correlative/predicted" 是准确的。

### 2.2 signature 基因的四条标准

原文（Methods p.e11）：

1. **ANOVA F 检验**：FDR adjusted P < 0.05（在所有簇之间）；
2. **Tukey HSD**：任一簇对 |log2 fold change| > 1 且 adjusted P < 0.01；
3. **AUC > 0.65**：该基因区分该簇与其余细胞的 ROC AUC；
4. **典型 marker**：文献中已知的细胞类型 marker。

四条缺一不可。第 4 条是**人工先验**，防止纯统计选出的基因没有生物学可解释性；复现时必须把 marker 列表落盘（`results/celltype_signatures.tsv` 的 `is_signature` 列 + `typical_marker` 列）。

### 2.3 相对丰度（signature score）

对 bulk 样本 $s$、细胞类型 $c$、signature 基因集合 $G_c$：

$$
z_{g,s} = \frac{\log_2(\mathrm{TPM}_{g,s}+1) - \mu_g}{\sigma_g}
$$

$$
A_{c,s} = \frac{1}{|G_c|}\sum_{g \in G_c} z_{g,s}
$$

其中 $\mu_g$、$\sigma_g$ 是该基因在**同一 cohort 内**所有样本的均值与标准差。z-score 的目的是去掉基因表达量级差异，使不同 signature 的丰度可比；但 z-score 也意味着**丰度是相对的**，不能解释为绝对细胞数。

### 2.4 基因–丰度相关

对每个基因 $g$ 与每个细胞类型 $c$，在 $N$ 个 bulk 样本上计算 Pearson 相关：

$$
r_{g,c} = \frac{\sum_{s=1}^{N}(x_{g,s}-\bar{x}_g)(A_{c,s}-\bar{A}_c)}
{\sqrt{\sum_s (x_{g,s}-\bar{x}_g)^2}\sqrt{\sum_s (A_{c,s}-\bar{A}_c)^2}}
$$

$r_{g,c}$ 高表示：当 $c$ 的丰度高时，$g$ 的表达也高。$g$ 可能是 $c$ 自己表达的（self-expressed），也可能是与 $c$ 共变的其它细胞表达的。

### 2.5 self-expressed 过滤

原文定义（Methods p.e11）：在 Smart-seq2 中，若基因 $g$ 在细胞类型 $c$ 中

1. 平均表达 > 1（log2 TPM），且
2. 表达细胞比例 > 20%，

则 $g$ 是 $c$ 的 self-expressed gene，把 $r_{g,c}$ 置 0。目的：避免"用 c 自己的 marker 证明 c 与 c 相关"的循环论证。

### 2.6 top 13 与富集分数

对每个细胞类型 $c$，在调整后的相关矩阵中取 $r'_{g,c}$ 最高的 13 个非自表达基因，记为 $T_c$。然后回到 Smart-seq2：

$$
z_{g,k} = \frac{\bar{e}_{g,k} - \mu_g^{sc}}{\sigma_g^{sc}}
$$

$$
ES_{k,c} = \frac{1}{13}\sum_{g \in T_c} z_{g,k}
$$

其中 $\bar{e}_{g,k}$ 是基因 $g$ 在细胞类型 $k$ 中的平均表达，$\mu_g^{sc}$、$\sigma_g^{sc}$ 是该基因在所有细胞类型中的均值与标准差。$ES_{k,c}$ 高表示：与 $c$ 丰度共变的 13 个基因，主要在 $k$ 中表达。再把 $ES_{k,c}$ 在 $k$ 上 z-score 化：

$$
Z_{k,c} = \frac{ES_{k,c} - \overline{ES_{\cdot,c}}}{\mathrm{sd}(ES_{\cdot,c})}
$$

$Z_{k,c} > 1.96$（单侧 95%）→ 连边 $c \to k$。原文用 1.96，正好对应正态单侧 97.5% 分位（双侧 95% 的一半）；复现时按原文，不改为 1.645 或 2.0。

### 2.7 为什么是 24 + 3

- 40 个免疫亚群中，只有 **24 个在 Smart-seq2 与 10x 两平台都组织富集（Ro/e>1）**；双平台一致可降低平台特异假象。
- 再加 3 个**肿瘤富集基质簇**：`hF01_Myofib-ACTA2`、`hF02_CAF-FAP`、`hE01_Endothelium-ACKR1`。原文明确说这三个是 tumor-enriched。
- 正常黏膜网络只用 normal-enriched 簇 + GTEx colon；肿瘤网络只用 tumor-enriched 簇 + TCGA COAD/READ。**不能把 normal 的簇丢进 tumor 网络**，否则丰度估计没有对应的 bulk 队列。

### 2.8 混杂与负对照（必须在报告中写）

| 混杂 | 机制 | 缓解 |
|---|---|---|
| **肿瘤纯度** | TCGA 样本的肿瘤细胞比例不同；纯度低 → 免疫 signature 高 → 假边 | 把纯度（ESTIMATE/ABSOLUTE）作协变量；或对丰度做纯度校正后敏感性分析 |
| **signature 共享** | 两个细胞类型共享 marker（如 `CD68` 在多个巨噬簇） | self-expressed 过滤 + 每簇独立 signature；报告共享基因数 |
| **样本构成** | TCGA COAD 与 READ 的 stage/部位不同 | 按 cohort 分层；Cox/相关中校正 stage |
| **生态相关** | 两个细胞类型被同一炎症状态驱动 | 负对照：随机 signature（打乱基因标签）应产生少量边；正常组织网络作对照 |
| **z-score 的相对性** | 不同 cohort 的 z-score 不可直接比较 | 肿瘤与正常分别 z-score，不合并 |

> ⚠️ **导师批注**：这一节是全文最容易被过度解读的地方。原文说 "predicted cell-cell interaction network"，不是 "interaction network"。复现报告的措辞必须是 **correlative / co-occurrence / predicted**，不能写成 "SPP1⁺ TAM 与 CAF 互作" 而定格为因果。

---

## ③ 原文逐句拆解

| 原文句（Methods / Results） | 中文转述 | 我们的操作 | 原文是否给参数 | 缺失项裁决 |
|---|---|---|---|---|
| "For 40 immune cell subtypes identified in human CRC, we assigned them into tumor or normal-enriched group based on their tissue enrichment scores quantified by Ro/e values (Ro/e > 1)" | 40 个免疫亚群按 Ro/e>1 分肿瘤/正常富集 | 从 `results/ro_e.tsv` 读；阈值 >1 | 给了阈值 | 无 |
| "Among them, 24 tissue-enriched clusters identified by both Smart-seq2 and 10x Genomics platforms were used" | 24 个双平台一致的组织富集簇 | 两平台 Ro/e 交集 | 给了 24，未列名单 | Table S4A 为准；否则按规则选（D11） |
| "we also included three tumor-enriched stromal cell subsets, including fibroblasts (hF01_Myofib-ACTA2 and hF02_CAF-FAP), and endothelial cells (hE01_Endothelium-ACKR1)" | 加 3 个基质簇 | 固定这三个标签 | 给了标签 | 无 |
| "The specific gene signatures for each subtype were derived from Smart-seq2 dataset based on ... FDR adjusted P-value of F test < 0.05" | ANOVA FDR<0.05 | 逐基因 `aov` + BH | 给了 | 无 |
| "any cluster pair showed significant difference in the HSD test (absolute log2 fold change > 1 and P-value adjusted by Tukey's HSD < 0.01)" | 任一簇对 Tukey | `TukeyHSD` 逐对 | 给了 | 无 |
| "AUC > 0.65" | AUC>0.65 | Wilcoxon AUC | 给了 | 无 |
| "typical marker genes described in literatures" | 典型 marker | 人工 marker 表 | 未列名单 | 由 04 号 marker 表 + 文献整理（D11） |
| "For tumor-enriched clusters, the TCGA COAD and READ data were used ... while the GTEx colon tissue data were used for normal-enriched clusters" | 肿瘤用 TCGA，正常用 GTEx | 分开计算 | 给了 | 无 |
| "estimated the relative abundance ... by the average expression of z-score normalized log-transformed expression of the cell type specific genes" | z-score 均值 | 公式 2.3 | 给了 | 无 |
| "we computed the Pearson correlation coefficient between the expression of each gene and the relative abundance of each cell subtype" | Pearson | 公式 2.4 | 给了 | 无 |
| "self-expressed genes ... average expression > 1; cell frequency of expression > 20%" | 自表达过滤 | 公式 2.5 | 给了 | 无 |
| "the top 13 highly correlated non-self-expressed genes were selected" | top13 | 每簇 13 个 | 给了 13 | 无 |
| "mean value of z-score transformed expression of these genes in each cell subtype as the enrichment score" | 富集分数 | 公式 2.6 | 给了 | 无 |
| "z-score transformed enrichment score > 1.96" | 连边阈值 | z>1.96 | 给了 | 无 |
| "correlative networks in tumor and normal mucosa were built separately" | 分开建网 | 两个网络 | 给了 | 无 |

---

## ④ 逐步可运行代码

> 运行环境：R 4.5.1 / Bioconductor 3.22 / Seurat 5.1.0 / glmnet 4.1-8。
> 所有脚本首行 `set.seed(20261018)`；输出目录不存在时 `dir.create(..., recursive=TRUE)`；缺输入时 `stop()`，不自动猜路径。

### 4.1 `code/06_01_signatures.R`：从 Smart-seq2 推导 signature

```r
#!/usr/bin/env Rscript
# code/06_01_signatures.R
# 用法: Rscript code/06_01_signatures.R <project_dir>
# 输入: processed/human_smartseq2_norm.rds, processed/human_smartseq2_annot.rds
# 输出: results/celltype_signatures.tsv, processed/signature_inputs.rds

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pROC)
})

set.seed(20261018)
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
stopifnot(dir.exists(file.path(proj, "processed")))
dir.create(file.path(proj, "results"), showWarnings = FALSE, recursive = TRUE)

sce <- readRDS(file.path(proj, "processed", "human_smartseq2_norm.rds"))
ann <- readRDS(file.path(proj, "processed", "human_smartseq2_annot.rds"))
common <- intersect(colnames(sce), colnames(ann))
if (length(common) < ncol(sce)) warning("Smart-seq2 norm/annot 细胞不完全一致，取交集：", length(common))
sce <- sce[, common]; ann <- ann[, common]
stopifnot(identical(colnames(sce), colnames(ann)))
colData(sce)$cluster_id <- ann$cluster_id
colData(sce)$tissue     <- ann$tissue

# 只用免疫 + 三个基质簇；排除恶性/上皮（避免肿瘤纯度主导）
keep_clusters <- c(
  grep("^hM", unique(sce$cluster_id), value = TRUE),
  grep("^hT", unique(sce$cluster_id), value = TRUE),
  grep("^hB", unique(sce$cluster_id), value = TRUE),
  grep("^hI", unique(sce$cluster_id), value = TRUE),
  "hF01_Myofib-ACTA2", "hF02_CAF-FAP", "hE01_Endothelium-ACKR1"
)
sce <- sce[, sce$cluster_id %in% keep_clusters]
sce <- sce[rowSums(assay(sce, "logtpm") > 0) >= 10, ]

expr <- as.matrix(assay(sce, "logtpm"))
clusters <- levels(factor(sce$cluster_id))
# 用 C1..Ck 作为 aov/Tukey 的安全标签，避免簇名中的连字符与对比名冲突
safe_levels <- paste0("C", seq_along(clusters))
cl <- factor(safe_levels[match(as.character(sce$cluster_id), clusters)],
             levels = safe_levels)

# ---- 1) 逐基因 ANOVA + Tukey（只保留在至少两组有表达的基因） ----
# 为控制耗时：先过滤在所有簇中平均表达 < 0.5 的基因
mean_by_cluster <- sapply(seq_along(clusters), function(i) {
  rowMeans(expr[, cl == safe_levels[i], drop = FALSE])
})
colnames(mean_by_cluster) <- clusters
keep_gene <- rowMaxs(mean_by_cluster) >= 0.5
expr <- expr[keep_gene, , drop = FALSE]
mean_by_cluster <- mean_by_cluster[keep_gene, , drop = FALSE]

anova_tukey <- function(g) {
  y <- expr[g, ]
  fit <- aov(y ~ cl)
  p_f <- summary(fit)[[1]][["Pr(>F)"]][1]
  tk  <- TukeyHSD(fit)$cl
  # 每个簇与其他簇的最大 |log2FC| 与最小 adj p
  # 安全标签 C1..Ck 使对比名形如 "C1-C2"，不会与簇名中的连字符冲突
  best <- sapply(seq_along(clusters), function(i) {
    k <- safe_levels[i]
    rn <- rownames(tk)
    hit <- grepl(paste0("(^|-)", k, "(-|$)"), rn)
    if (!any(hit)) return(c(max_abs = 0, min_p = 1))
    d <- tk[hit, "diff"]; p <- tk[hit, "p adj"]
    c(max_abs = max(abs(d)), min_p = min(p))
  })
  data.frame(
    gene = g,
    p_anova = p_f,
    max_abs_log2FC = best["max_abs", ],
    min_tukey_p = best["min_p", ],
    cluster_id = clusters,
    row.names = NULL
  )
}

# 逐基因循环；15k 基因 × 40 簇在 8 核上约 10–20 min
res_list <- lapply(rownames(expr), anova_tukey)
res <- bind_rows(res_list)
# 每个簇内跨基因 BH（原文的 FDR 是 gene-level；按簇分别校正更接近其口径）
res <- res %>%
  group_by(cluster_id) %>%
  mutate(FDR_anova = p.adjust(p_anova, method = "BH")) %>%
  ungroup()

# ---- 2) AUC：每个簇 vs 其余 ----
auc_one <- function(g, k) {
  y <- expr[g, ]
  pos <- y[cl == k]; neg <- y[cl != k]
  if (length(pos) < 10 || length(neg) < 10) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(response = c(rep(1, length(pos)), rep(0, length(neg))),
                                 predictor = c(pos, neg), quiet = TRUE)))
}

# 只对通过 ANOVA/Tukey 的 gene×cluster 计算 AUC，节省时间
cand <- res %>% filter(FDR_anova < 0.05, max_abs_log2FC > 1, min_tukey_p < 0.01)
cand$auc <- mapply(function(g, k) auc_one(g, k), cand$gene, cand$cluster_id)

# ---- 3) 典型 marker（人工整理；完整表见 04 号 §2.2） ----
marker_panel <- list(
  "hM01_Mast-TPSAB1" = c("TPSAB1","TPSB2","CPA3","KIT","MS4A2"),
  "hM02_pDC-LILRA4" = c("LILRA4","IL3RA","GZMB","TCF4"),
  "hM03_cDC2-CD1C" = c("CD1C","FCER1A","CLEC10A","CD1E"),
  "hM04_cDC1-BATF3" = c("BATF3","XCR1","CLEC9A","IRF8"),
  "hM05_Mono-CD14" = c("CD14","S100A8","S100A9","VCAN","LYZ"),
  "hM06_Mono-CD16" = c("FCGR3A","CX3CR1","MS4A7","LST1"),
  "hM07_Mono-CD14CD16" = c("CD14","FCGR3A","S100A8","S100A9"),
  "hM08_Macro-NLRP3" = c("NLRP3","IL1B","EREG","CSF1R"),
  "hM09_Macro-PLTP" = c("PLTP","LYVE1","IL10","FOLR2"),
  "hM10_Macro-IL1B" = c("IL1B","C1QA","C1QB","APOE"),
  "hM11_Monolike-FCN1" = c("FCN1","CD14","S100A8","C1QC","SPP1"),
  "hM12_TAM-C1QC" = c("C1QA","C1QB","C1QC","TREM2","MERTK","CD80"),
  "hM13_TAM-SPP1" = c("SPP1","MARCO","VEGFA","FN1","SDC2"),
  "hF01_Myofib-ACTA2" = c("ACTA2","TAGLN","MYL9","SERPING1"),
  "hF02_CAF-FAP" = c("FAP","COL1A1","COL3A1","DCN","MMP2"),
  "hE01_Endothelium-ACKR1" = c("ACKR1","CDH5","ERG","VWF","PLVAP"),
  "hT01_CD4-CCR7" = c("CCR7","SELL","LEF1","TCF7"),
  "hT02_CD4-ANXA1" = c("ANXA1","IL7R","GPR183"),
  "hT03_CD4-GNLY" = c("GNLY","GZMA","GZMB"),
  "hT04_CD4-TCF7" = c("TCF7","LEF1","SELL"),
  "hT05_CD4-CXCR6" = c("CXCR6","ITGAE","ZNF683"),
  "hT06_CD4-CXCR5" = c("CXCR5","BCL6","PDCD1"),
  "hT07_CD4-GZMK" = c("GZMK","GZMA","EOMES"),
  "hT08_CD4-IL23R" = c("IL23R","RORC","CCR6"),
  "hT09_CD4-CXCL13" = c("CXCL13","PDCD1","CTLA4"),
  "hT10_CD4-FOXP3" = c("FOXP3","CTLA4","IL2RA"),
  "hT11_CD4-CTLA4" = c("CTLA4","FOXP3","TNFRSF4"),
  "hT12_CD8-LEF1" = c("LEF1","SELL","TCF7"),
  "hT13_CD8-GPR183" = c("GPR183","CCR7","IL7R"),
  "hT14_CD8-CX3CR1" = c("CX3CR1","FCGR3A","FGFBP2"),
  "hT15_CD8-GZMK" = c("GZMK","GZMA","EOMES"),
  "hT16_CD8-CD6" = c("CD6","TRAT1","CD3D"),
  "hT17_CD8-CD160" = c("CD160","KIR2DL4","ZNF683"),
  "hT18_CD8-LAYN" = c("LAYN","HAVCR2","PDCD1","TOX"),
  "hB01_PlasmaB-IgG" = c("MZB1","JCHAIN","IGHG1"),
  "hB02_GALTB-IgA" = c("GALTB","IGHA1","IGHA2"),
  "hB03_FollicularB-IgD" = c("IGHD","MS4A1","CD19"),
  "hB04_FollicularB-MS4A1" = c("MS4A1","CD79A","CD37"),
  "hB05_GCBCell-LRMP" = c("LRMP","TCL1A","BANK1"),
  "hI01_NK-CD16" = c("FCGR3A","NKG7","GNLY"),
  "hI02_NK-GZMK" = c("GZMK","NKG7","EOMES"),
  "hI03_NK-CD103" = c("ITGAE","ZNF683","NKG7"),
  "hI04_ILC3-SLC4A10" = c("SLC4A10","KIT","IL23R","RORC")
)

cand$typical_marker <- mapply(function(g, k) {
  m <- marker_panel[[k]]
  if (is.null(m)) FALSE else g %in% m
}, cand$gene, cand$cluster_id)

cand <- cand %>%
  mutate(is_signature = FDR_anova < 0.05 &
           max_abs_log2FC > 1 &
           min_tukey_p < 0.01 &
           !is.na(auc) & auc > 0.65 &
           typical_marker)

write.table(cand, file.path(proj, "results", "celltype_signatures.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(expr = expr, cl = cl, cl_orig = factor(sce$cluster_id),
             mean_by_cluster = mean_by_cluster),
        file.path(proj, "processed", "signature_inputs.rds"))

cat("signature gene×cluster 条目:", nrow(cand),
    "；is_signature=TRUE:", sum(cand$is_signature), "\n")
```

> **实现说明**：`marker_panel` 的键必须与 03 号的 `cluster_id` **完全一致**（`hM01_Mast-TPSAB1` 这种"下划线分大类、连字符分亚型"的写法）；不要用 `gsub()` 批量替换，否则会把 `hM01_` 也换成 `hM01-`，导致 `marker_panel[[k]]` 取不到值。典型 marker 的完整表以 04 号 `results/myeloid_markers.tsv` 与 Table S2 为准。

### 4.2 `code/06_02_abundance.R`：TCGA/GTEx 的 signature 丰度

```r
#!/usr/bin/env Rscript
# code/06_02_abundance.R
# 用法: Rscript code/06_02_abundance.R <project_dir>
# 输入: raw/tcga_coadread/expr.tsv, raw/tcga_coadread/clinical.tsv
#       raw/gtex_colon/expr.tsv, raw/gtex_colon/sample_anno.tsv
#       results/celltype_signatures.tsv
# 输出: results/abundance_tcga.tsv, results/abundance_gtex.tsv

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
})
set.seed(20261018)
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
dir.create(file.path(proj, "results"), showWarnings = FALSE, recursive = TRUE)

read_expr <- function(path) {
  # 期望: 第一列 gene symbol，其余列 sample；值可以是 log2(TPM+1) 或 raw TPM
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  rownames(dt) <- dt[[1]]; dt <- dt[, -1, drop = FALSE]
  m <- as.matrix(dt)
  if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)  # 粗略判断 raw TPM
  m
}

expr_tcga <- read_expr(file.path(proj, "raw", "tcga_coadread", "expr.tsv"))
expr_gtex <- read_expr(file.path(proj, "raw", "gtex_colon", "expr.tsv"))
sig <- read.delim(file.path(proj, "results", "celltype_signatures.tsv"))
sig <- sig[sig$is_signature, ]

# 只保留肿瘤富集/正常富集簇（由 03 号 Ro/e 表给出；这里从 signature 表读 cluster_id）
# 约定：abundance 对所有 signature 簇都算；下游 06_04 再用 Ro/e 选边
compute_abundance <- function(expr, sig, cohort) {
  genes <- intersect(rownames(expr), unique(sig$gene))
  expr <- expr[genes, , drop = FALSE]
  # 每个基因在所有样本内 z-score；sd=0 时置 0
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  out <- sapply(unique(sig$cluster_id), function(k) {
    g <- intersect(sig$gene[sig$cluster_id == k], rownames(z))
    if (length(g) < 3) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE])
  })
  out <- as.data.frame(out, check.names = FALSE)
  out$sample <- colnames(expr)
  out$cohort <- cohort
  out
}

ab_tcga <- compute_abundance(expr_tcga, sig, "TCGA")
ab_gtex <- compute_abundance(expr_gtex, sig, "GTEx")

write.table(ab_tcga, file.path(proj, "results", "abundance_tcga.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ab_gtex, file.path(proj, "results", "abundance_gtex.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("TCGA samples:", nrow(ab_tcga), "; GTEx samples:", nrow(ab_gtex), "\n")
```

> **注意**：TCGA 与 GTEx 必须**分别** z-score，不能合并；否则 cohort 批次会主导丰度。TCGA 的 `expr.tsv` 建议用 Xena 的 `TOIL RSEM norm_count` 或 `TPM`，在 `01` 号记录版本与下载日期。若使用 `log2(norm_count+1)`，公式中的 TPM 替换为 norm_count，但 z-score 后的相对丰度解释不变；契约中必须写明尺度（`PROJECT_SPEC.md` §4）。

### 4.3 `code/06_03_correlation.R`：基因–丰度相关矩阵

```r
#!/usr/bin/env Rscript
# code/06_03_correlation.R
# 用法: Rscript code/06_03_correlation.R <project_dir>
# 输入: raw/tcga_coadread/expr.tsv, raw/gtex_colon/expr.tsv,
#       results/abundance_tcga.tsv, results/abundance_gtex.tsv,
#       results/celltype_signatures.tsv, processed/signature_inputs.rds
# 输出: results/correlation_gene_cluster.tsv, processed/correlation_matrices.rds

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
})
set.seed(20261018)
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
stopifnot(file.exists(file.path(proj, "results", "celltype_signatures.tsv")))

read_expr <- function(path) {
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  rownames(dt) <- dt[[1]]; dt <- dt[, -1, drop = FALSE]
  m <- as.matrix(dt)
  if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)
  m
}

expr <- list(
  Tumor  = read_expr(file.path(proj, "raw", "tcga_coadread", "expr.tsv")),
  Normal = read_expr(file.path(proj, "raw", "gtex_colon", "expr.tsv"))
)
ab <- list(
  Tumor  = read.delim(file.path(proj, "results", "abundance_tcga.tsv"), check.names = FALSE),
  Normal = read.delim(file.path(proj, "results", "abundance_gtex.tsv"), check.names = FALSE)
)
sig <- read.delim(file.path(proj, "results", "celltype_signatures.tsv"))
sig <- sig[sig$is_signature, ]
inputs <- readRDS(file.path(proj, "processed", "signature_inputs.rds"))

# self-expressed: Smart-seq2 平均表达 >1 且表达比例 >20%
self_expr <- function(g, k) {
  if (!g %in% rownames(inputs$expr)) return(FALSE)
  cells <- which(inputs$cl_orig == k)
  if (length(cells) < 10) return(FALSE)
  v <- inputs$expr[g, cells]
  isTRUE(mean(v, na.rm = TRUE) > 1 && mean(v > 0, na.rm = TRUE) > 0.20)
}
self_tab <- expand.grid(gene = unique(sig$gene), cluster_id = unique(sig$cluster_id),
                        stringsAsFactors = FALSE)
self_tab$self <- mapply(self_expr, self_tab$gene, self_tab$cluster_id)

cor_one <- function(tissue) {
  e <- expr[[tissue]]; a <- ab[[tissue]]
  stopifnot(all(a$sample %in% colnames(e)))
  e <- e[, a$sample, drop = FALSE]
  # 只保留在至少 20% 样本中表达的基因，降低方差为 0 的噪声
  e <- e[rowMeans(e > 0) >= 0.2, , drop = FALSE]
  clusters <- setdiff(colnames(a), c("sample", "cohort"))
  # Pearson: 基因 × 样本 vs 丰度 × 样本
  rmat <- cor(t(e), as.matrix(a[, clusters, drop = FALSE]), method = "pearson")
  pmat <- matrix(NA_real_, nrow(rmat), ncol(rmat), dimnames = dimnames(rmat))
  n <- ncol(e)
  tstat <- rmat * sqrt((n - 2) / pmax(1 - rmat^2, 1e-12))
  pmat <- 2 * pt(-abs(tstat), df = n - 2)
  # self-expressed 置 0
  idx <- self_tab[self_tab$self, ]
  for (i in seq_len(nrow(idx))) {
    g <- idx$gene[i]; k <- idx$cluster_id[i]
    if (g %in% rownames(rmat) && k %in% colnames(rmat)) rmat[g, k] <- 0
  }
  list(r = rmat, p = pmat, n = n, tissue = tissue)
}

res <- lapply(c("Tumor", "Normal"), cor_one)
names(res) <- c("Tumor", "Normal")
saveRDS(res, file.path(proj, "processed", "correlation_matrices.rds"))

long <- bind_rows(lapply(res, function(x) {
  as.data.frame(as.table(x$r)) %>%
    rename(gene = Var1, cluster_id = Var2, r = Freq) %>%
    mutate(p = as.vector(x$p), n = x$n, tissue = x$tissue)
}))
# 契约列 is_self：与 self_tab 对齐（self_tab$self 为 TRUE 表示自表达）
long <- long %>%
  left_join(self_tab %>% select(gene, cluster_id, is_self = self),
            by = c("gene", "cluster_id"))
write.table(long, file.path(proj, "results", "correlation_gene_cluster.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("correlation rows:", nrow(long), "\n")
```

### 4.4 `code/06_04_network.R`：top13 → 富集分数 → 网络边

```r
#!/usr/bin/env Rscript
# code/06_04_network.R
# 用法: Rscript code/06_04_network.R <project_dir>
# 输入: processed/correlation_matrices.rds, processed/signature_inputs.rds,
#       results/ro_e.tsv, results/celltype_signatures.tsv
# 输出: results/correlation_network_edges.tsv, results/network_nodes.tsv

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(matrixStats)
})
set.seed(20261018)
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
stopifnot(file.exists(file.path(proj, "processed", "correlation_matrices.rds")))

cors <- readRDS(file.path(proj, "processed", "correlation_matrices.rds"))
inputs <- readRDS(file.path(proj, "processed", "signature_inputs.rds"))
roe <- read.delim(file.path(proj, "results", "ro_e.tsv"))
sig <- read.delim(file.path(proj, "results", "celltype_signatures.tsv"))
sig <- sig[sig$is_signature, ]

# 组织特异的富集簇：肿瘤网络用 Tumor 富集簇 + 3 个基质簇；正常网络用 Normal 富集簇
stromal <- c("hF01_Myofib-ACTA2", "hF02_CAF-FAP", "hE01_Endothelium-ACKR1")
cluster_sets <- lapply(c(Tumor = "Tumor", Normal = "Normal"), function(tt) {
  cl <- roe %>% filter(tissue == tt, ro_e > 1) %>% pull(cluster_id) %>% unique()
  if ("platform" %in% colnames(roe)) {
    ok <- roe %>% filter(cluster_id %in% cl, tissue == tt) %>%
      group_by(cluster_id) %>%
      summarise(n_platform = n_distinct(platform), .groups = "drop") %>%
      filter(n_platform >= 2) %>% pull(cluster_id)
    cl <- intersect(cl, ok)
  }
  cl
})
cluster_sets$Tumor <- unique(c(cluster_sets$Tumor, stromal))
stopifnot(length(cluster_sets$Tumor) > 0, length(cluster_sets$Normal) > 0)

# Smart-seq2 每个簇的平均表达（基因 × 簇）
mean_expr <- inputs$mean_by_cluster
# 每个基因在所有簇中的 z-score
z_expr <- t(scale(t(mean_expr)))
z_expr[!is.finite(z_expr)] <- 0

build_network <- function(tissue) {
  x <- cors[[tissue]]
  r <- x$r
  # 只保留该组织的富集簇（且在相关矩阵中）
  clusters <- intersect(cluster_sets[[tissue]], colnames(r))
  edges <- list()
  for (c in clusters) {
    # top13 非自表达基因（r 已被置 0；按 r 降序）
    ord <- order(r[, c], decreasing = TRUE)
    top <- rownames(r)[ord][seq_len(min(13, nrow(r)))]
    top <- intersect(top, rownames(z_expr))
    if (length(top) < 5) next
    # ES: 每个簇 k 的 top 基因 z-score 均值
    es <- colMeans(z_expr[top, , drop = FALSE])
    # z-score ES across clusters
    z_es <- (es - mean(es)) / sd(es)
    hit <- names(z_es)[z_es > 1.96]
    hit <- setdiff(hit, c)  # 去掉自己
    hit <- intersect(hit, cluster_sets[[tissue]])  # 只保留纳入网络的节点
    if (length(hit) == 0) next
    edges[[length(edges) + 1]] <- data.frame(
      source_cluster = c, target_cluster = hit,
      n_top_genes = length(top), enrichment_z = as.numeric(z_es[hit]),
      edge = TRUE, tissue = tissue, stringsAsFactors = FALSE
    )
  }
  bind_rows(edges)
}

edges <- bind_rows(lapply(c("Tumor", "Normal"), build_network))
if (nrow(edges) == 0) stop("网络无任何边：检查 ro_e.tsv、signature 与相关矩阵是否对齐")
# 补上未连边的节点（便于画图）
all_nodes <- bind_rows(lapply(names(cluster_sets), function(tt) {
  expand.grid(source_cluster = cluster_sets[[tt]], target_cluster = cluster_sets[[tt]],
              tissue = tt, stringsAsFactors = FALSE)
}))
edges_full <- all_nodes %>%
  left_join(edges, by = c("source_cluster", "target_cluster", "tissue")) %>%
  mutate(edge = ifelse(is.na(edge), FALSE, edge),
         n_top_genes = ifelse(is.na(n_top_genes), 0L, n_top_genes))

write.table(edges_full, file.path(proj, "results", "correlation_network_edges.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 无向网络：每个边对 source/target 各计一次，避免只数出边导致 Hub 被低估
nodes <- bind_rows(
  edges_full %>% filter(edge) %>% transmute(cluster_id = source_cluster, tissue),
  edges_full %>% filter(edge) %>% transmute(cluster_id = target_cluster, tissue)
) %>%
  count(cluster_id, tissue, name = "degree") %>%
  arrange(tissue, desc(degree))
write.table(nodes, file.path(proj, "results", "network_nodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("edges (TRUE):", sum(edges_full$edge), "\n")
```

### 4.5 `code/06_05_plot_network.R`：Fig. 3E / Fig. S5C

```r
#!/usr/bin/env Rscript
# code/06_05_plot_network.R
# 用法: Rscript code/06_05_plot_network.R <project_dir>
# 输入: results/correlation_network_edges.tsv
# 输出: figures/Fig3E_network.pdf/.png, figures/FigS5C_network_normal.pdf/.png

suppressPackageStartupMessages({
  library(dplyr); library(igraph); library(ggraph); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
proj <- ifelse(length(args) >= 1, args[1], ".")
dir.create(file.path(proj, "figures"), showWarnings = FALSE)

edges <- read.delim(file.path(proj, "results", "correlation_network_edges.tsv"))
nodes <- read.delim(file.path(proj, "results", "network_nodes.tsv"))
major_type <- function(x) {
  dplyr::case_when(
    grepl("^hM", x) ~ "Myeloid",
    grepl("^hT", x) ~ "T cell",
    grepl("^hB", x) ~ "B cell",
    grepl("^hI", x) ~ "ILC/NK",
    grepl("^hF", x) ~ "Fibroblast",
    grepl("^hE", x) ~ "Endothelium",
    TRUE ~ "Other"
  )
}

plot_one <- function(tissue, out_pdf, out_png) {
  e <- edges %>% filter(edge, tissue == !!tissue)
  n <- nodes %>% filter(tissue == !!tissue)
  g <- graph_from_data_frame(e[, c("source_cluster","target_cluster","enrichment_z")],
                             vertices = n, directed = FALSE)
  V(g)$type <- major_type(V(g)$name)
  V(g)$size <- 3 + 0.6 * degree(g)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = enrichment_z), alpha = 0.35, color = "grey40") +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    geom_node_point(aes(size = size, color = type)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.4) +
    scale_size_identity() +
    theme_void() +
    labs(title = paste0(tissue, " correlative network"),
         color = "Cell type") +
    theme(legend.position = "right")
  ggsave(out_pdf, p, width = 9, height = 7)
  ggsave(out_png, p, width = 9, height = 7, dpi = 300)
}

plot_one("Tumor",
         file.path(proj, "figures", "Fig3E_network.pdf"),
         file.path(proj, "figures", "Fig3E_network.png"))
plot_one("Normal",
         file.path(proj, "figures", "FigS5C_network_normal.pdf"),
         file.path(proj, "figures", "FigS5C_network_normal.png"))
```

### 4.6 手算示例（3 基因 × 3 样本）

假设某 signature 含 $g_1,g_3$，bulk 表达（log2 TPM）为：

| 基因 | s1 | s2 | s3 | 均值 | sd | z(s1) | z(s2) | z(s3) |
|---|---|---|---|---|---|---|---|---|
| g1 | 5 | 6 | 7 | 6 | 1 | −1 | 0 | 1 |
| g3 | 8 | 8 | 8 | 8 | 0 | 0 | 0 | 0 |

丰度：$A = \frac{1}{2}[(-1,0,1)+(0,0,0)] = (-0.5, 0, 0.5)$。

若 $g_2=(2,2,3)$，均值 2.333、sd 0.577；$A$ 均值 0、sd 0.5：

$$
r = \frac{\frac{1}{2}[(-0.5)(-0.333)+0\cdot(-0.333)+0.5\cdot0.667]}{0.5 \times 0.577}
= \frac{0.25}{0.2885} \approx 0.866
$$

若 $g_2$ 在 signature 中（self-expressed），则把 $r$ 置 0；否则它可能进入 top13，用于连边。这个例子说明：**z-score 的 sd=0 处理必须显式置 0，否则会产生 NaN**；代码里 `z[!is.finite(z)] <- 0` 就是为此。

---

## ⑤ 参数解释与调参空间

| 参数 | 原文取值 | 可选范围 | 改变后果 | 本项目裁决 |
|---|---|---|---|---|
| 纳入簇 | 24 免疫 + 3 基质 | 20–30 + 3 | 少 → 网络稀疏；多 → 假边增加 | 原文 24+3；Table S4A 为准（D11） |
| Ro/e 阈值 | >1 | 0.8–2 | 低 → 纳入弱富集；高 → 只留强富集 | 1.0 |
| signature FDR | <0.05 | 0.01–0.1 | 严 → 基因少、AUC 高；松 → 基因多、噪声 | 0.05 |
| Tukey | \|log2FC\|>1, adj p<0.01 | 0.5–2 / 0.001–0.05 | 同上 | 原文 |
| AUC | >0.65 | 0.6–0.8 | 严 → 特异但少；松 → 多但共享 | 0.65 |
| 典型 marker | 人工 | — | 决定第 4 条是否通过 | 04 号 marker 表 + 文献 |
| top N | 13 | 5–30 | 少 → 富集分数不稳；多 → 混入非特异 | 13（原文） |
| 连边 z | >1.96 | 1.645–2.58 | 低 → 边多；高 → 边少 | 1.96（原文） |
| 相关方法 | Pearson | Spearman | Pearson 对 outlier 敏感；Spearman 更稳健 | Pearson 主路径；Spearman 敏感性 |
| z-score 范围 | 每 cohort | 全样本/分 cohort | 合并会引入 cohort 批次 | 分 cohort |
| 纯度校正 | 原文未做 | ESTIMATE/ABSOLUTE | 纯度低 → 免疫丰度假高 | 主路径不做（照原文）；敏感性做 |

---

## ⑥ 预期结果与 QC 检查点

> 对应 `14` 号 §六 Q27–Q31。文献值来自 Fig. 3E、S5C、S8A；未报告具体边数/度数的地方标注「文献未报告」。

| # | 指标 | 文献值 | 允许范围 | Go/No-Go | 我方值（S6） |
|---|---|---|---|---|---|
| Q27 | 纳入簇数 | 24 + 3 | 27 ± 2 | ✅/⚠️ | |
| Q28 | 网络 Hub | TAM 与 cDC 度数最高 | 前 5 名含至少 2 个 TAM/cDC | ✅/⚠️/❌ | |
| Q29 | SPP1⁺ 边 | CAF、myofibroblast、endothelium | 至少 2 类 | ✅/⚠️/❌ | |
| Q30 | C1QC⁺ 边 | T 细胞亚群 | 至少 2 个 T 簇 | ✅/⚠️/❌ | |
| Q31 | 正常黏膜网络 | B–T–DC（淋巴滤泡） | 方向一致 | ✅/⚠️ | |
| Q31b | 肿瘤网络边数 | 文献未报告 | 报告实际值；与 S5C 视觉比较 | 记录 | |
| Q31c | 负对照（随机 signature） | 文献未报告 | 边数应显著少于真实网络 | ✅/⚠️ | |
| Q31d | 纯度敏感性 | 文献未报告 | 校正后 Hub 不变 | 记录 | |

**判读规则**：

- 如果 Hub 不是 TAM/cDC，先检查 signature 是否混入上皮/恶性基因；再检查 TCGA 纯度；最后检查 self-expressed 过滤是否生效。
- 如果 SPP1⁺ 没有 CAF/内皮边，检查 `hF01/hF02/hE01` 是否被纳入；检查 `SPP1`、`FN1`、`SDC2` 是否在 signature；检查 `r` 的符号。
- 如果正常网络没有 B–T–DC，检查 GTEx colon 是否包含淋巴滤泡样本；检查 `hB03/hB04` 与 Tfh（`hT06`）是否 normal-enriched。

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| `scale()` 产生 NaN | 某基因在所有样本中表达相同（sd=0） | `z[!is.finite(z)] <- 0`；记录被置 0 的基因数 |
| 相关矩阵全是 NA | 样本顺序不匹配；`a$sample` 与 `colnames(e)` 不一致 | 用 `match()` 显式对齐；`stopifnot(all(a$sample %in% colnames(e)))` |
| 网络边爆炸 | top13 取了绝对值最大而非正相关；或 z 阈值用了 1.0 | 只取 `r` 降序；z>1.96；检查 `order()` 的 `decreasing=TRUE` |
| 网络没有边 | signature 基因太少；top13 与 z_expr 交集 <5；ES 的 sd=0 | 打印每簇 top 基因数与 `length(top)`；检查 `z_expr` 是否全 0 |
| Hub 是上皮/恶性 | signature 未排除 `hC01–hC04`/`hE02–hE06` | `keep_clusters` 只保留免疫 + 3 基质；检查 `grep("^hM|^hT|^hB|^hI")` |
| 肿瘤/正常网络混在一起 | 没按 `tissue` 分组 | `build_network("Tumor")` 与 `"Normal"` 分开；丰度表也分开 |
| TCGA 与 GTEx 的 z-score 不可比 | 合并 z-score | 分 cohort z-score；报告中说明 |
| 相关被 outlier 驱动 | 单个高表达样本 | 报告 Spearman 敏感性；检查杠杆点（hat values） |
| 边方向无法解释 | 相关性无方向 | 只画无向边；不要画箭头；措辞用 "co-occurrence" |
| 复现的边与 Fig. 3E 不一致 | 原文用 Table S4A 的 signature；我们重新推导 | 优先用 Table S4A；若不可得，记录 D11 并做敏感性 |

---

## ⑧ 与文献图表的对应关系

| 原文图 | 内容 | 本模块产物 | 说明 |
|---|---|---|---|
| Fig. 3E | 肿瘤相关性网络 | `figures/Fig3E_network.pdf/.png` | 节点大小=度数；颜色=大类 |
| Fig. S5A | myofib/CAF–endothelium 的相关性示例 | `results/correlation_gene_cluster.tsv` 子集 | 可单独画散点 + 回归 |
| Fig. S5B | 高相关基因的富集热图 | 可由 `06_04` 的 `z_expr` 与 top13 列表派生（模块私有中间对象，不单独登记契约） | 行=簇，列=TCGA 样本/基因 |
| Fig. S5C | 正常黏膜网络 | `figures/FigS5C_network_normal.pdf/.png` | B–T–DC 结构 |
| Fig. S8A | 分析流程图 | 本文件 §2 | 五步流程 |
| Table S4A | signature 基因列表 | `results/celltype_signatures.tsv` | 若原文表可得，以其为准 |

---

## 附：本模块的"一句话答辩"

**我们不是在做细胞通讯，而是在做"细胞共变"：用 scRNA-seq 定义细胞指纹，用 bulk 估计指纹丰度，用相关找共变伙伴，再用 scRNA 把伙伴基因放回细胞类型。它能提出假设（SPP1⁺–CAF/内皮、C1QC⁺–T），但不能证明互作；证明需要空间/功能实验（07/13）。**
