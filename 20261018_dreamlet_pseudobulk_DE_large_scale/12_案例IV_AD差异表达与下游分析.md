# 12 案例 IV：AD 差异表达与下游分析

> **目标文献**：Hoffman GE, et al. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* Nature Communications (2026) 17:9597. DOI: 10.1038/s41467-026-75680-8
> **模块对应原文**：Fig. 5C–H、Fig. 6A–F；Supplementary Figs. 7–8；Results「Modeling multiple sources of expression variation in a large Alzheimer's Disease cohort」「Large effect up-regulation of PTPRG in microglia of Alzheimer's disease cases」；作者代码 `PsychAD_r0/PsychAD_r0_analysis.Rmd`
> **前置依赖**：`11_案例III_AD队列数据生成与上游处理.md`（`processed/sce_psychad_annot.rds`）或直接 `raw/psychad_ad/*.h5ad`（D1）
> **复现难度 / 预计工时 / 计算资源**：★★★★☆ / 5–8 人日 / 内存 64 GB；全流程 1–2 天（32 核）

---

## ① 科学问题与逻辑位置

- **上游**：PsychAD snRNA-seq（1.4 M nuclei / 586 samples / 299 donors / 22 细胞类型）。
- **下游**：Fig. 5–6 全部面板；13 号的跨方法/跨队列对照；14 号的 QC 门槛 G3。
- **本模块回答**：**"AD 在每个细胞类型里改变了哪些基因？效应多大？microglia 的 PTPRG 上调是不是全局最强信号？为什么需要 299 个供体才能看见？"**

这是全文的"旗舰案例"，也是三个案例中**唯一同时具备大样本、重复测量、高维批次**的设计，因此最能体现 dreamlet 的价值。

---

## ② 背景原理

### 2.1 AD 的细胞生物学

| 细胞类型 | 在 AD 中的角色 | 关键基因 |
|---|---|---|
| **小胶质细胞（Micro_PVM）** | 淀粉样斑块清除、神经炎症、tau 传播 | `PTPRG`, `DUSP10`, `APOE`, `TREM2`, `CSF1R` |
| **星形胶质细胞（Astro）** | 反应性增生、谷氨酸稳态 | `GFAP`, `AQP4`, `SLC1A2` |
| **少突胶质细胞（Oligo）** | 髓鞘形成/修复 | `PLP1`, `MBP`, `MAG` |
| **OPC** | 髓鞘再生 | `PDGFRA`, `VCAN` |
| **兴奋性神经元（EN）** | 突触丢失、兴奋毒性 | `SLC17A7`, `RBFOX3` |
| **抑制性神经元（IN）** | 网络平衡（PVALB/SST/VIP） | `GAD1`, `GAD2`, `PVALB`, `SST` |
| **血管/脑膜（Endo/PC/VLMC）** | 血脑屏障、AD 血管假说 | `CLDN5`, `PDGFRB`, `SLC14A1` |

### 2.2 关键发现背后的生物学

- **PTPRG**（protein tyrosine phosphatase receptor gamma）：受体型蛋白酪氨酸磷酸酶；原文报告 microglia 中 log2FC=1.59、p=2.94e-28，是"效应量与显著性都站得住"的少数基因之一。原文讨论：PTPRG 是炎症标志物，且是**三项前人 microglia 转录组研究中两项都上调**的三个基因之一（另两个是 APOE、DPYD）；AD GWAS 未在该位点找到风险变异 → 提示这是**反应性上调**，可能随疾病阶段变化。
- **PDE10A**：磷酸二酯酶，纹状体/皮层中等棘神经元高表达；AD 中在 13/14 个神经元亚型下调——神经元丢失或转录下调。
- **DUSP10**：双特异性磷酸酶，p38MAPK 通路的负反馈调控；microglia 中上调，与 p38MAPK cascade 激活一致。
- **p38MAPK**：应激/炎症通路；在 microglia 的炎症反应中起核心作用（原文引文 40）。

> ⚠️ **导师批注**：PTPRG 的"最大效应"是**统计上的**，不是"因果上的"。原文自己写："its role in the molecular etiology of AD is unclear"、"upregulation is reactive and may vary with the stage and progression"。复现时不要把"最大 logFC"写成"驱动基因"（L0 §12 相关性≠因果）。

### 2.3 研究设计：为什么需要 299 个供体

AD 的效应量普遍温和（原文：仅 101 个基因被 AD 解释 >5% 方差）。在 22 个细胞类型 × ~15k 基因的多重检验下，要达到 5% study-wide FDR：

- 小样本（<50）时，即使真实效应存在，检验也会被"方差估计噪声 + 多重校正"淹没；
- 大样本时，EB 收缩的 $d_0$ 更稳定、每基因的 $d_g$ 更大，moderated t 的零分布更准；
- 原文 Fig. 5G–H 正是展示：**每个受试者的核数越多 → 技术重复一致性越高 → DEG 越多**。

### 2.4 模型与公式（作者代码）

```r
# 权重公式（processAssays）
~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + (1|Dx)

# DE 公式（dreamlet）
~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + Dx

# VP 公式（fitVarPart）
~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + (1|Dx)
```

| 项 | 类型 | 含义 |
|---|---|---|
| `SubID` | 随机 | 受试者（299） |
| `poolID` | 随机 | 10x pool（数量以实际列为准，D18） |
| `Sex` | 随机 | 性别（在 DE 里作随机效应，作者代码如此） |
| `scale(Age)` | 固定 | 年龄（标准化） |
| `Dx` | 固定（DE）/ 随机（VP） | AD vs Control |

> ⚠️ **导师批注**：把 `Sex` 放进随机效应是作者的建模选择（D5），不是唯一正确做法。性别只有 2 个水平，理论上更适合固定效应；用随机效应会让它与其他变量共享收缩。复现时必须照做，但敏感性分析要跑"Sex 固定"版本，并在方法学说明。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "the dataset comprised 586 samples from 299 donors over 60 years of age and 1.4 M nuclei. Of these donors, 150 had Alzheimer's disease (AD), and 149 were age-matched neurotypical controls." | 队列规模 | 数据核对 | 给了数字 | D12 |
| "Cell cluster annotation identified 22 cell types, including eight subtypes of excitatory neurons and 6 subtypes of inhibitory neurons." | 22 类 | 注释 | 给了数字 | — |
| "Analysis of pseudobulk counts for microglia (Micro_PVM) shows the strong mean-variance trend." | voom | 03 号 | — | — |
| "the median fraction of variation explained by subject is 38.4%" | Subject 方差 | 06 号 | 给了数字 | — |
| "only 101 genes explaining more than 5% of the variance" | AD 效应温和 | 06 号 | 给了数字 | — |
| "932 genes have >5% variance explained across the pools" | 批次效应 | 06 号 | 给了数字 | — |
| "20.2% of the variation in PTPRG is explained by Alzheimer's status while 97.3% of variance in NXPE1 is explained by sample pool." | 代表基因 | 06 号 | 给了数字 | — |
| "microglia have 1037 differentially expressed genes passing the study-wide 5% FDR cutoff." | 1037 DEG | 05 号 | 给了数字 | D6 |
| "PTPRG stands out with a log2 fold change of 1.59 and p-value of 2.94e-28." | PTPRG | 05 号 | 给了数字 | — |
| "PDE10A is significantly down-regulated in AD in 13 of 14 neuron subtypes." | PDE10A | 跨 assay | 给了数字 | — |
| "microglia have a specific up-regulation of the p38MAPK cascade … This pathway includes DUSP10, which is one of the top upregulated genes in microglia." | p38MAPK | 07 号 | 给了方向 | D9 |
| "increasing the average number of nuclei per subject is directly related to an increase in the technical reproducibility" | Fig. 5G | VP vs 核数 | 未给统计量 | 用 Subject 方差分数 |
| "the number of nuclei per subject impacts the read count per subject, the number of genes that pass a minimum expression cutoff, and the number of genes that are found to be differentially expressed" | Fig. 5H | 三联动 | 未给回归 | 线性回归 + p |

---

## ④ 逐步可运行代码

### 4.1 pseudobulk + 权重 + DE + VP（AD 主流程）

```r
# code/12_01_dreamlet_ad.R
suppressPackageStartupMessages({
  library(dreamlet); library(SingleCellExperiment); library(BiocParallel)
})

sce <- readRDS("processed/sce_psychad_annot.rds")   # 或从 raw h5ad 读
sce$Dx <- factor(sce$Dx_AD, levels = c("Control","AD"))

pb <- aggregateToPseudoBulk(sce, assay = "counts",
                            cluster_id = "subclass", sample_id = "Channel",
                            BPPARAM = SnowParam(8))
saveRDS(list(pb = pb, cellCounts = cellCounts(pb)), "processed/pb_psychad.rds")

# 第一级 Poisson 权重
w0 <- lapply(assayNames(pb), function(k) {
  y <- as.matrix(assay(pb, k)); w <- y + 0.5; w / rowMeans(w)
})
names(w0) <- assayNames(pb)
saveRDS(w0, "processed/weights_poisson_psychad.rds")

form_w <- ~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + (1|Dx)
proc <- processAssays(pb, form_w, min.cells = 5, min.count = 10,
                      weightsList = w0, BPPARAM = SnowParam(12))
saveRDS(proc, "processed/proc_psychad.rds")

fit <- dreamlet(proc, ~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + Dx,
                BPPARAM = SnowParam(12))
saveRDS(fit, "processed/fit_psychad.rds")
writeLines(coefNames(fit), "logs/12_coefnames_psychad.txt")

vp <- fitVarPart(proc, form_w, BPPARAM = SnowParam(12))
saveRDS(vp, "processed/vp_psychad.rds")
```

### 4.2 核心数字：1037 DEG、PTPRG 1.59

```r
# code/12_02_key_numbers.R
suppressPackageStartupMessages({ library(dreamlet); library(dplyr) })

fit <- readRDS("processed/fit_psychad.rds")
stopifnot("DxAD" %in% coefNames(fit))

tt <- bind_rows(lapply(assayNames(fit), function(k) {
  x <- as.data.frame(topTable(fit[[k]], coef = "DxAD", number = Inf, sort.by = "none"))
  if (nrow(x) == 0) return(NULL)
  x$ID <- rownames(x); x$assay <- k; x$coef <- "DxAD"
  x
}))
# ★ study-wide BH 必须**按 coef 分组**（D6），不能把多个 coef 混在一起校正
tt <- tt |> dplyr::group_by(coef) |>
  dplyr::mutate(adj.P.Val.studywide = p.adjust(P.Value, method = "BH")) |>
  dplyr::ungroup()
write.table(tt, "results/topTable_psychad_DxAD.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

mg <- tt[tt$assay == "Micro_PVM", ]
cat("Micro_PVM genes:", nrow(mg), " DEG (study-wide):",
    sum(mg$adj.P.Val.studywide < 0.05), "\n")

ptprg <- mg[mg$ID == "PTPRG", ]
print(ptprg[, c("ID","logFC","t","P.Value","adj.P.Val","adj.P.Val.studywide")])
cat("expected: logFC≈1.59, P≈2.94e-28\n")

# PDE10A 在多少神经元亚型中显著下调
neuron <- tt[grepl("^(EN|IN)_", tt$assay), ]
pde <- neuron[neuron$ID == "PDE10A", ]
cat("PDE10A significant (study-wide) in", sum(pde$adj.P.Val.studywide < 0.05 & pde$logFC < 0),
    "of", nrow(pde), "neuron subtypes\n")
```

### 4.3 Fig. 5G–H：核数与效力

```r
# code/12_03_nuclei_power.R
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

pb   <- readRDS("processed/pb_psychad.rds")
cc   <- pb$cellCounts
vp   <- as.data.frame(readRDS("processed/vp_psychad.rds"))
tt   <- read.delim("results/topTable_psychad_DxAD.tsv")

# 每个细胞类型：**平均每受试者**核数（cc 的行是 Channel，不是 SubID！）
# 先把 Channel 的细胞数按 SubID 求和，再对受试者取均值
sample_subject <- setNames(pb$pb$SubID, colnames(pb$pb))
nuc_per_subj <- sapply(assayNames(pb$pb), function(k) {
  per <- tapply(cc[, k], sample_subject[rownames(cc)], sum, na.rm = TRUE)
  mean(per, na.rm = TRUE)
})

# Subject 方差分数（技术重复一致性）：作者代码用该细胞类型内 Subject 的均值
subj_var <- vp |>
  filter(!is.na(SubID)) |>
  group_by(assay) |>
  summarize(mean_subject = mean(SubID), median_subject = median(SubID),
            .groups = "drop") |>
  mutate(mean_nuclei = nuc_per_subj[assay])

# DEG 比例
deg <- tt |>
  group_by(assay) |>
  summarize(n_genes = n(), frac_deg = mean(adj.P.Val.studywide < 0.05), .groups = "drop")
df <- left_join(subj_var, deg, by = "assay") |>
  mutate(mean_nuclei = as.numeric(mean_nuclei))

# 回归
fit1 <- lm(mean_subject ~ log10(mean_nuclei), data = df)   # 也可用 median_subject 做敏感性
fit2 <- lm(frac_deg ~ log10(mean_nuclei), data = df)
cat("Fig5G: subject variance ~ log10(nuclei): p =", summary(fit1)$coefficients[2, 4], "\n")
cat("Fig5H: frac DEG ~ log10(nuclei): p =", summary(fit2)$coefficients[2, 4], "\n")

p1 <- ggplot(df, aes(mean_nuclei, mean_subject, size = n_genes)) +
  geom_point(alpha = 0.7) + scale_x_log10() + theme_classic() +
  labs(x = "Mean nuclei per subject", y = "Variance explained by subject (mean)")
ggsave("figures/Fig5G_nuclei_reproducibility.png", p1, width = 5, height = 4, dpi = 300)

p2 <- ggplot(df, aes(mean_nuclei, frac_deg)) +
  geom_point() + scale_x_log10() + theme_classic() +
  labs(x = "Mean nuclei per subject", y = "Fraction DEG (study-wide FDR<5%)")
ggsave("figures/Fig5H_nuclei_deg.png", p2, width = 5, height = 4, dpi = 300)
```

### 4.4 Fig. 6 出图

```r
# code/12_04_figures_ad.R
suppressPackageStartupMessages({ library(dreamlet); library(ggplot2) })

fit <- readRDS("processed/fit_psychad.rds")
ctorder <- c("EN_L2_3_IT","EN_L3_5_IT_1","EN_L3_5_IT_2","EN_L3_5_IT_3","EN_L5_6_NP",
             "EN_L6_CT","EN_L6_IT","EN_NF","IN_ADARB2","IN_LAMP5","IN_PVALB",
             "IN_PVALB_CHC","IN_SST","IN_VIP","Oligo","OPC","Astro","Micro_PVM",
             "CD8_T","PC","VLMC","Endo")

# Fig. 6A：microglia 火山图
p1 <- plotVolcano(fit, coef = "DxAD", assays = "Micro_PVM")
ggsave("figures/Fig6A_volcano_Micro_PVM.pdf", p1, width = 5, height = 5)
ggsave("figures/Fig6A_volcano_Micro_PVM.png", p1, width = 5, height = 5, dpi = 300)

# Fig. 6B：跨细胞类型 z 统计量热图
genes <- c("PTPRG","DUSP10","NPNT","DPYD","VGF","SPRY4-AS1","PDE10A","NCAM2","PCDH7")
p2 <- plotGeneHeatmap(fit, coef = "DxAD", genes = genes, assays = ctorder)
ggsave("figures/Fig6B_gene_heatmap.pdf", p2, width = 10, height = 4)
ggsave("figures/Fig6B_gene_heatmap.png", p2, width = 10, height = 4, dpi = 300)

# Fig. 6C：PTPRG 森林图
p3 <- plotForest(fit, coef = "DxAD", gene = "PTPRG", assays = ctorder)
ggsave("figures/Fig6C_forest_PTPRG.pdf", p3, width = 6, height = 8)

# Fig. 6D：PDE10A 森林图
p4 <- plotForest(fit, coef = "DxAD", gene = "PDE10A", assays = ctorder)
ggsave("figures/Fig6D_forest_PDE10A.pdf", p4, width = 6, height = 8)
```

### 4.5 基因集分析（Fig. 6E–F）

```r
# code/12_05_zenith_ad.R
suppressPackageStartupMessages({ library(dreamlet); library(zenith); library(GSEABase); library(msigdbr) })

fit <- readRDS("processed/fit_psychad.rds")
# ★ 作者代码用 GO（不是 Hallmark）：p38 = GO1900744，synapse assembly = GO 相关条目
go.gs <- get_GeneOntology("BP", to = "SYMBOL")
go.gs <- go.gs[sapply(go.gs, function(x) length(geneIds(x))) < 1000]
stopifnot("GO1900744" %in% names(geneIds(go.gs)))

res <- zenith_gsa(fit, go.gs, coefs = "DxAD", progressbar = FALSE)
write.table(res, "results/zenith_psychad_DxAD_GO_BP.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 敏感性分析：Hallmark（D9 的第二套）
hall <- msigdbr(species = "Homo sapiens", collection = "H")
hall.lst <- split(hall$gene_symbol, hall$gs_name)
hall.gsc <- GeneSetCollection(mapply(function(g, nm) {
  GeneSet(unique(g), setName = nm, geneIdType = SymbolIdentifier())
}, hall.lst, names(hall.lst), SIMPLIFY = FALSE))
res_h <- zenith_gsa(fit, hall.gsc, coefs = "DxAD", progressbar = FALSE)
write.table(res_h, "results/zenith_psychad_DxAD_Hallmark.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

p <- plotZenithResults(res, 5, 3)
ggsave("figures/Fig6E_zenith_AD.pdf", p, width = 14, height = 10)
ggsave("figures/Fig6E_zenith_AD.png", p, width = 14, height = 10, dpi = 300)

# Fig. 6F：p38MAPK 通路基因热图
# ⚠️ Hallmark 没有 p38MAPK 集合；用 GO BP（GOBP_P38MAPK_CASCADE）
go.gs <- get_GeneOntology("BP", to = "SYMBOL")
stopifnot("GOBP_P38MAPK_CASCADE" %in% names(geneIds(go.gs)))
p38_genes <- geneIds(go.gs[["GOBP_P38MAPK_CASCADE"]])
p38_genes <- intersect(p38_genes, rownames(topTable(fit[["Micro_PVM"]], coef = "DxAD", number = Inf)))
p2 <- plotGeneHeatmap(fit, coef = "DxAD", genes = p38_genes, assays = "Micro_PVM")
ggsave("figures/Fig6F_p38MAPK_genes.png", p2, width = 8, height = 6, dpi = 300)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| DE 公式 | §2.4 | Sex 固定、加 `(1|round_num)` | 影响方差分配与 p 值（D5） |
| `min.cells` / `min.count` | 5 / 10 | 3/5、10/10 | 细胞类型数与基因数 |
| `subclass` vs `subtype` | `subclass`（22 类） | `subtype`（更细） | 更细 → 每类细胞少、效力低 |
| FDR 口径 | study-wide | 细胞类型内 | DEG 数差一个数量级（D6） |
| `Sex` 的角色 | 随机（作者） | 固定 | 敏感性分析 |
| 基因集 | Hallmark | GO | 通路结论（D9） |
| 核数-效力回归 | log10 线性 | 样条 | 只影响图，不影响结论 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| 样本 / 供体 / 细胞 | 586 / 299 / 1.4 M | ±5% | ❌ D12 |
| 细胞类型数 | 22 | 可映射 | ❌ |
| microglia DEG（study-wide） | 1037 | ±10% | ❌ |
| PTPRG log2FC | 1.59 | ±10% | ❌ |
| PTPRG p | 2.94e-28 | 同数量级 | ⚠️ |
| PTPRG 的 AD 方差分数 | 20.2% | ±3 pp | ⚠️ |
| NXPE1 的 pool 方差分数 | 97.3% | ±3 pp | ⚠️ |
| Subject 中位方差 | 38.4% | ±3 pp | ❌ |
| AD >5% 的基因数 | 101 | ±15% | ⚠️ |
| pool >5% 的基因数 | 932 | ±15% | ⚠️ |
| PDE10A 显著下调的神经元亚型 | 13/14 | ±1 | ⚠️ |
| DUSP10 在 microglia 上调 | 是 | — | ❌ |
| p38MAPK 通路显著 | 是 | — | ⚠️ |
| 核数 vs Subject 方差 | 正相关、p<0.05 | — | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| PTPRG 不显著 | 样本版本不同 / FDR 口径不同 | 核对 D12/D6；先看 logFC 方向 |
| 1037 变 3000+ | 用了细胞类型内 BH | 改 study-wide |
| 22 类对不上 | 上游注释不同 | 映射到 22 类；不要硬调 |
| `DxAD` 系数不存在 | `Dx` 因子水平不是 Control/AD | `levels(sce$Dx)` |
| VP 的 Sex 分数异常高 | 性别与 AD 分组不平衡 | 报告性别分布；敏感性分析 |
| 核数-效力回归不显著 | 细胞类型太少/噪声大 | 报告 p 与趋势；不要过度解释 |
| `plotGeneHeatmap` 灰块多 | 基因在某细胞类型未过表达过滤 | 正常；图注说明 |
| 结果与作者差 2%–10% | 数据 release/QC 阈值不同 | 记录在 gap 表 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 5C（voom） | `figures/Fig5C_voom_Micro_PVM.*` | 03 号生成 |
| Fig. 5D（VP） | `figures/Fig5D_vp_AD.*` | 06 号生成 |
| Fig. 5E（PTPRG/NXPE1） | `figures/Fig5E_vp_examples.*` | 06 号 |
| Fig. 5F（GC 相关） | `figures/Fig5F_gc_pool_correlation.png` | 06 号 |
| Fig. 5G（核数-可重复性） | `figures/Fig5G_nuclei_reproducibility.png` | §4.3 |
| Fig. 5H（核数-DEG） | `figures/Fig5H_nuclei_deg.png` | §4.3 |
| Fig. 6A（火山） | `figures/Fig6A_volcano_Micro_PVM.*` | §4.4 |
| Fig. 6B（跨细胞类型热图） | `figures/Fig6B_gene_heatmap.*` | §4.4 |
| Fig. 6C（PTPRG 森林） | `figures/Fig6C_forest_PTPRG.*` | §4.4 |
| Fig. 6D（PDE10A 森林） | `figures/Fig6D_forest_PDE10A.*` | §4.4 |
| Fig. 6E（GSA） | `figures/Fig6E_zenith_AD.*` | §4.5 |
| Fig. 6F（p38MAPK） | `figures/Fig6F_p38MAPK_genes.png` | §4.5 |
| Fig. S7–S8 | `results/12_ad_summary.csv` | 汇总 |

> **交付物**：`processed/pb_psychad.rds`、`proc_psychad.rds`、`fit_psychad.rds`、`vp_psychad.rds`、`results/topTable_psychad_DxAD.tsv`、`results/zenith_psychad_DxAD_Hallmark.tsv`、`figures/Fig5C–H`、`figures/Fig6A–F`。完成后进入 `13_方法对照与跨队列一致性.md`。
