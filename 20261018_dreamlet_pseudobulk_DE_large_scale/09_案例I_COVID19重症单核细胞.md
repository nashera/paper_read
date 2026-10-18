# 09 案例 I：COVID-19 重症单核细胞

> **目标文献**：Hoffman GE, et al. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* Nature Communications (2026) 17:9597. DOI: 10.1038/s41467-026-75680-8
> **模块对应原文**：Fig. 3A–H；Supplementary Fig. 4；Results「Dreamlet uncovers cell-type-specific response associated with COVID-19 severity」；作者代码 `COVID_combat/COVID_combat_major.Rmd`、`COVID_combat_minor.Rmd`
> **前置依赖**：`02`–`07` 的管线；`raw/combat_covid19/` 的 h5ad 与临床变量表
> **复现难度 / 预计工时 / 计算资源**：★★★★☆ / 3–5 人日 / 内存 64 GB；全队列 DE 1–3 h（16 核）

---

## ① 科学问题与逻辑位置

- **上游**：COMBAT 血液 scRNA-seq（674 K cells / 110 donors）+ 临床严重度。
- **下游**：Fig. 3 全部面板；13 号的跨方法对照；14 号的 QC 门槛 G2/G3。
- **本模块回答**：**"COVID-19 的疾病严重度在单核细胞里留下怎样的细胞类型特异转录印记？它与脓毒症（Sepsis）有何不同？"**

选择单核细胞的原因（原文）：它们是炎症反应的核心，且在 COVID-19 中表现出最显著的亚群重塑（经典单核细胞 `cMono`、非经典 `ncMono`、循环型 `cMono.cyc`、血小板相关 `cMono.PLT`、树突状 `cDC`；以作者 `major_subset` 标签为准）。

---

## ② 背景原理

### 2.1 生物学背景

| 细胞类型 | 标志基因 | 功能 | 在 COVID 中的变化 |
|---|---|---|---|
| **经典单核细胞（cMono / CD14⁺）** | `CD14`, `LYZ`, `S100A8/9` | 吞噬、炎症因子分泌 | 重症时扩增、HLA-DR 下降、干扰素反应增强 |
| **非经典单核细胞（ncMono / CD16⁺）** | `FCGR3A`, `MS4A7`, `CDKN1C` | 血管巡逻、组织修复 | 重症时扩增；原文 Fig. 3A 的重点之一 |
| **循环型（cMono.cyc）** | 细胞周期基因（`MKI67`, `TOP2A`） | 增殖态单核 | 数量少、变异大 |
| **血小板相关（cMono.PLT）** | `PPBP`, `PF4` | 与血小板黏附的单核 | 技术/生物学混杂 |
| **树突状细胞（DC）** | `FCER1A`, `CLEC10A` | 抗原提呈 | 重症时减少 |
| **血小板（Platelet）** | `PPBP`, `PF4` | 凝血 | COVID 血栓并发症相关 |

原文的关键基因：

- **CLU（clusterin）**：补体抑制、脂质运输、应激反应；在 COVID-19 与脓毒症中均升高（Fig. 3F）。
- **NFKBIA（IκBα）**：NF-κB 的负反馈抑制因子；在 COVID-19 中上调、脓毒症中不上调——提示 COVID 特异的 NF-κB 通路激活（Fig. 3F）。
- **TNFα/NF-κB 通路**：COVID 各严重度均激活，脓毒症不激活（Fig. 3G–H）。

### 2.2 队列设计：为什么必须处理"一个 donor 多个样本"

COMBAT 的设计是**纵向/多部位采样**：同一患者可能有多个 `scRNASeq_sample_ID`（不同时间点/不同中心）。作者代码做了两步：

1. **每位 donor 只保留最重症的一个样本**（`which.max(Source)`）；
2. **只保留 Oxford 中心**的样本（因为 `Source` 与 `Institute` 混杂：Flu 只在 St George's 收集）。

> ⚠️ **导师批注**：这两步是**结论口径**层面的选择（D10）。若保留同一 donor 的多个样本，就必须把 `donor_id` 放进随机效应（`(1|donor_id)`），否则 FPR 会膨胀。作者选择"每人一个样本"来简化设计——这不是唯一正确做法，但必须与原文一致。**不要**在不声明的情况下改成"全样本 + donor 随机效应"。

### 2.3 `Source` 的六个水平与参考组

作者代码的 `coef_array`：

```r
c("SourceCOVID_HCW_MILD", "SourceCOVID_MILD", "SourceCOVID_SEV",
  "SourceCOVID_CRIT", "SourceSepsis")
```

对应六个水平（参考组为 `Healthy`）：

| 水平 | 含义 |
|---|---|
| `Healthy` | 健康对照（参考） |
| `COVID_HCW_MILD` | 医护人员轻症（职业暴露） |
| `COVID_MILD` | 轻症（非住院） |
| `COVID_SEV` | 重症（住院） |
| `COVID_CRIT` | 危重症（ICU） |
| `Sepsis` | 脓毒症（疾病对照） |

> ⚠️ **先核对文件里的实际顺序**：`levels(sce$Source)` 若已是有序因子，沿用文件顺序（这才是作者 `which.max(Source)` 的依据），并在 D10 记录；若为字符型，再按本节显式设定。**不要**把本节顺序当成原文事实。
>
> ⚠️ **必须显式设定因子水平顺序**。若 `Source` 是字符型，R 会按字母序排参考组（`COVID_CRIT` 可能成为参考），所有系数名与符号都会变。作者代码用 `which.max(Source)` 选最重症，也依赖因子的**有序编码**——复现时要先 `factor(Source, levels = c("Healthy","COVID_HCW_MILD","COVID_MILD","COVID_SEV","COVID_CRIT","Sepsis"))`，再确认 `which.max` 的语义。

### 2.4 模型

- **权重公式**：`~ Age + (1|sex) + (1|Source)`（`Source` 作随机效应，吸收六种状态的方差）。
- **DE 公式**：`~ Source + Age + sex`（`Source` 作固定效应，逐水平检验）。
- **VP 公式**：同权重公式，输出每个基因被 `Source`/`Age`/`sex` 解释的比例。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "The COMBAT Consortium collected blood from patients hospitalized with COVID-19 or sepsis, as well as from non-hospitalized COVID-19 patients and healthy controls." | 四类人群 | 数据背景 | — | — |
| "We performed dreamlet analysis on single-cell transcriptome data of 674K cells from 110 donors." | 674K/110 | 样本筛选后 | 给了数字 | D10 |
| "we focus on monocyte subtypes and related cell types due to their key role in inflammatory response." | 聚焦单核 | `major_subset` | — | — |
| "Analysis of pseudobulk counts for classical monocytes shows a strong mean-variance trend in the log2 counts which dreamlet models using two-step precision weights." | voom 趋势 | 03 号 | — | — |
| "Variance partitioning analysis … median of 11.4% and 3250 of the 12,472 genes … >25%." | VP | 06 号 | 给了数字 | — |
| "Comparing gene expression in classical monocytes between patients with mild COVID-19 and healthy controls identified 1824 differentially expressed genes." | 1824 DEG | `coef="SourceCOVID_MILD"` | 给了数字 | D6 |
| "CLU expression increasing with COVID-19 severity and sepsis compared to healthy controls while NFKBIA upregulated in COVID-19 but not sepsis patients." | 两基因 | 森林图 | 给了方向 | — |
| "TNFα signaling via NF-κB showed activation in all COVID-19 disease states but not in patients with sepsis." | 通路 | zenith | 给了方向 | D9 |

---

## ④ 逐步可运行代码

### 4.1 读入、清洗、合并临床变量

```r
# code/09_01_prepare_covid.R
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(dplyr)
})

sce <- readH5AD("raw/combat_covid19/covid_counts.h5ad", use_hdf5 = TRUE)
counts(sce) <- assay(sce, "X")          # 作者做法：原始计数在 X
rownames(sce) <- rowData(sce)$feature_name   # 用 SYMBOL

# 去双胞、去未分类
sce <- sce[, sce$GEX_region != "E: Doublets"]
sce <- sce[, sce$major_subset != "nan"]
sce$major_subset <- droplevels(factor(sce$major_subset))

# 临床变量
clin <- read.delim("raw/combat_covid19/COMBAT_CLINVAR_for_processed.txt")
clin <- clin[clin$scRNASeq_sample_ID %in% sce$scRNASeq_sample_ID, ]
idx  <- match(sce$scRNASeq_sample_ID, clin$scRNASeq_sample_ID)
colData(sce) <- cbind(colData(sce), clin[idx, ])

# 显式因子水平（见 §2.3）
sce$Source <- factor(sce$Source,
  levels = c("Healthy","COVID_HCW_MILD","COVID_MILD","COVID_SEV","COVID_CRIT","Sepsis"))
stopifnot(!any(is.na(sce$Source)))
```

### 4.2 每位 donor 一个样本 + 仅 Oxford

```r
# code/09_02_sample_selection.R
cd <- as.data.frame(colData(sce))

# 每位 donor 取 Source 编码最大的样本（最重症）
sel <- cd |>
  distinct(donor_id, scRNASeq_sample_ID, Source) |>
  group_by(donor_id) |>
  slice_max(order_by = as.integer(Source), n = 1, with_ties = FALSE) |>
  ungroup()

sce_sub <- sce[, sce$scRNASeq_sample_ID %in% sel$scRNASeq_sample_ID]
sce_sub <- sce_sub[, sce_sub$Institute == "Oxford"]
sce_sub$scRNASeq_sample_ID <- droplevels(factor(sce_sub$scRNASeq_sample_ID))

cat("donors:", length(unique(sce_sub$donor_id)),
    " samples:", length(unique(sce_sub$scRNASeq_sample_ID)),
    " cells:", ncol(sce_sub), "\n")
saveRDS(sce_sub, "processed/sce_covid_sub.rds")
```

### 4.3 pseudobulk + 权重 + DE + VP

```r
# code/09_03_dreamlet_covid.R
suppressPackageStartupMessages({
  library(dreamlet); library(SingleCellExperiment); library(BiocParallel)
})

sce <- readRDS("processed/sce_covid_sub.rds")

pb <- aggregateToPseudoBulk(sce, assay = "counts",
                            cluster_id = "major_subset",
                            sample_id  = "scRNASeq_sample_ID",
                            BPPARAM = SnowParam(8))
saveRDS(list(pb = pb, cellCounts = cellCounts(pb)), "processed/pb_covid.rds")

# 第一级 Poisson 权重
w0 <- lapply(assayNames(pb), function(k) {
  y <- as.matrix(assay(pb, k)); w <- y + 0.5
  w / rowMeans(w)
})
names(w0) <- assayNames(pb)
saveRDS(w0, "processed/weights_poisson_covid.rds")

form_w <- ~ Age + (1|sex) + (1|Source)
proc <- processAssays(pb, form_w, min.cells = 3, min.count = 3,
                      weightsList = w0, BPPARAM = SnowParam(8))
saveRDS(proc, "processed/proc_covid.rds")

fit <- dreamlet(proc, ~ Source + Age + sex, BPPARAM = SnowParam(8))
saveRDS(fit, "processed/fit_covid.rds")
writeLines(coefNames(fit), "logs/09_coefnames_covid.txt")

vp <- fitVarPart(proc, form_w, BPPARAM = SnowParam(8))
saveRDS(vp, "processed/vp_covid.rds")
```

### 4.4 结果提取与关键数字

```r
# code/09_04_results_covid.R
suppressPackageStartupMessages({ library(dreamlet); library(dplyr); library(qvalue) })

fit <- readRDS("processed/fit_covid.rds")
coefs <- c("SourceCOVID_HCW_MILD","SourceCOVID_MILD","SourceCOVID_SEV",
           "SourceCOVID_CRIT","SourceSepsis")

tt <- bind_rows(lapply(coefs, function(cf) {
  bind_rows(lapply(assayNames(fit), function(k) {
    x <- as.data.frame(topTable(fit[[k]], coef = cf, number = Inf, sort.by = "none"))
    if (nrow(x) == 0) return(NULL)
    x$ID <- rownames(x); x$assay <- k; x$coef <- cf; x$SE <- x$logFC / x$t
    x
  }))
}))
# ★ study-wide BH 必须**按 coef 分组**（D6），不能把多个 coef 混在一起校正
tt <- tt |> dplyr::group_by(coef) |>
  dplyr::mutate(adj.P.Val.studywide = p.adjust(P.Value, method = "BH")) |>
  dplyr::ungroup()
write.table(tt, "results/topTable_covid_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 关键数字
cmono <- tt[tt$assay == "cMono", ]
cat("cMono expressed genes:", nrow(cmono[cmono$coef == "SourceCOVID_MILD", ]), "\n")
cat("mild vs healthy DEG (study-wide):",
    sum(cmono$adj.P.Val.studywide[cmono$coef == "SourceCOVID_MILD"] < 0.05), "\n")

# VP 中位与 >25%
vp <- as.data.frame(readRDS("processed/vp_covid.rds"))
x <- vp[vp$assay == "cMono", ]
cat("median Source variance:", median(x$Source), "\n")
cat("genes >25% Source variance:", sum(x$Source > 0.25), "of", nrow(x), "\n")
```

### 4.5 Fig. 3 出图

```r
# code/09_05_figures_covid.R
suppressPackageStartupMessages({ library(dreamlet); library(ggplot2) })

fit  <- readRDS("processed/fit_covid.rds")
proc <- readRDS("processed/proc_covid.rds")

# Fig. 3B：voom 趋势
p <- plotVoom(proc, assays = "cMono")
ggsave("figures/Fig3B_voom_cMono.pdf", p, width = 5, height = 5)
ggsave("figures/Fig3B_voom_cMono.png", p, width = 5, height = 5, dpi = 300)

# Fig. 3C：VP
p2 <- plotVarPart(sortCols(readRDS("processed/vp_covid.rds")), assays = "cMono")
ggsave("figures/Fig3C_vp_cMono.pdf", p2, width = 5, height = 4)
ggsave("figures/Fig3C_vp_cMono.png", p2, width = 5, height = 4, dpi = 300)

# Fig. 3E：火山图
p3 <- plotVolcano(fit, coef = "SourceCOVID_MILD", assays = "cMono")
ggsave("figures/Fig3E_volcano_mild.pdf", p3, width = 5, height = 5)
ggsave("figures/Fig3E_volcano_mild.png", p3, width = 5, height = 5, dpi = 300)

# Fig. 3F：CLU / NFKBIA 森林图
p4 <- plotForest(fit, coef = "SourceCOVID_MILD", gene = "CLU", assays = "cMono")
ggsave("figures/Fig3F_forest_CLU.png", p4, width = 5, height = 4, dpi = 300)
```

### 4.6 基因集分析

```r
# code/09_06_zenith_covid.R
suppressPackageStartupMessages({ library(dreamlet); library(zenith); library(GSEABase) })

fit <- readRDS("processed/fit_covid.rds")
go.gs <- get_GeneOntology(to = "SYMBOL")
go.gs <- go.gs[sapply(go.gs, function(x) length(geneIds(x))) < 1000]

coefs <- c("SourceCOVID_HCW_MILD","SourceCOVID_MILD","SourceCOVID_SEV",
           "SourceCOVID_CRIT","SourceSepsis")
res <- zenith_gsa(fit, go.gs, coefs = coefs, progressbar = FALSE)
write.table(res, "results/zenith_covid_GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

p <- plotZenithResults(res, 5, 3)
ggsave("figures/Fig3G_gsa_covid.pdf", p, width = 14, height = 16)
ggsave("figures/Fig3G_gsa_covid.png", p, width = 14, height = 16, dpi = 300)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| 样本选择 | 每 donor 最重症 + Oxford | 全样本 + donor 随机效应 | 改变样本量与混杂结构（D10） |
| `min.cells` / `min.count` | 3 / 3 | 5/5 | 影响保留的细胞类型与基因数 |
| `Source` 水平顺序 | §2.3 | 字母序 | 参考组变 → 系数名/符号变 |
| `sex` 在 DE 中的角色 | 固定效应 | 随机效应 | 原文 DE 用固定；VP 用随机 |
| 基因集 | GO（<1000 基因） | Hallmark | 影响通路结论（D9） |
| FDR 口径 | study-wide | 细胞类型内 | DEG 数差一个数量级（D6） |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| donor 数 / 细胞数 | 110 / 674 K | ±10% | ❌ 查样本筛选（D10） |
| `cMono` 表达基因数 | 12,472 | ±10% | ❌ |
| `cMono` 疾病中位方差 | 11.4% | ±1.5 pp | ❌ |
| >25% 的基因数 | 3250 | ±10% | ❌ |
| mild vs healthy DEG | 1824 | ±10% | ❌ |
| CLU：随严重度升高 | 方向一致 | — | ❌ |
| NFKBIA：COVID 上调、Sepsis 不上调 | 方向一致 | — | ❌ |
| TNFα/NF-κB：COVID 激活、Sepsis 不激活 | 方向一致 | — | ❌ |
| 单核细胞亚群数 | 与原文一致 | — | ⚠️ |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 系数名不是 `SourceCOVID_MILD` | 参考组不是 `Healthy` | 显式 `factor(levels=...)` |
| `which.max(Source)` 选错样本 | `Source` 是字符/未排序因子 | 先转有序因子，或手工映射严重度 0–5 |
| 细胞数远少于 674K | 漏了 `major_subset != "nan"` 或双胞过滤 | 按 4.1 逐步做 |
| `Institute` 列不存在 | 临床表列名不同 | `colnames(clin)` 检查；用 `grep("institute", ..., ignore.case=TRUE)` |
| 某单核亚群缺失 | `major_subset` 标签不同（`cMono` vs `Classical`） | `table(sce$major_subset)`；按作者标签映射 |
| VP 的 `Source` 分数为 0 | `Source` 在 processAssays 里被当固定效应 | 检查权重公式 |
| 结果与原文差很多 | 用了全部 Institute | 只保留 Oxford（D10） |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 对应代码 |
|---|---|---|
| Fig. 3A（单核亚群 UMAP） | `figures/Fig3A_mono_umap.*` | `plotProjection(sce, "X_umap", "major_subset")` |
| Fig. 3B（voom 趋势） | `figures/Fig3B_voom_cMono.*` | `plotVoom()` |
| Fig. 3C（VP 小提琴） | `figures/Fig3C_vp_cMono.*` | `plotVarPart()` |
| Fig. 3D（代表基因） | `figures/Fig3D_vp_examples.*` | `plotPercentBars()` |
| Fig. 3E（火山） | `figures/Fig3E_volcano_mild.*` | `plotVolcano()` |
| Fig. 3F（森林） | `figures/Fig3F_forest_*` | `plotForest()` |
| Fig. 3G（通路热图） | `figures/Fig3G_gsa_covid.*` | `plotZenithResults()` |
| Fig. 3H（TNFα 基因） | `figures/Fig3H_tnfa_genes.*` | `plotGeneHeatmap()` |
| Fig. S4（全细胞类型） | `results/09_de_summary_covid.csv` | 汇总 |

> **交付物**：`processed/pb_covid.rds`、`proc_covid.rds`、`fit_covid.rds`、`vp_covid.rds`、`results/topTable_covid_all.tsv`、`results/zenith_covid_GO.tsv`、`figures/Fig3*`。完成后进入 `10_案例II_前列腺癌骨转移.md`。
