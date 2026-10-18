# 07 基因集分析：camera 与 zenith

> **目标文献**：Hoffman GE, et al. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* Nature Communications (2026) 17:9597. DOI: 10.1038/s41467-026-75680-8
> **模块对应原文**：Fig. 3G–H（COVID 通路热图、TNFα/NF-κB 基因热图）、Fig. 4G–H（MHC-II、蛋白折叠 chaperone）、Fig. 6E–F（AD 通路、p38MAPK/DUSP10）；Methods「Gene set analysis using the full spectrum of test statistics」（引文 33 camera、74 MSigDB）
> **前置依赖**：`05_差异表达分析流程与结果提取.md`（`processed/fit_<dataset>.rds`）、`06_方差分解与批次效应量化.md`（可选，用于解释）
> **复现难度 / 预计工时 / 计算资源**：★★★☆☆ / 1–2 人日 / 内存 16 GB；AD 全队列 0.5–1 h

---

## ① 科学问题与逻辑位置

- **上游**：`dreamletResult`（每个基因的 moderated t）+ 基因集（GO / MSigDB）。
- **下游**：`results/zenith_<dataset>_<geneset>.tsv`、`figures/Fig*G/H_*`。
- **本模块回答**：**"哪些通路/基因集在疾病中整体上调或下调？"** ——并且要避免两个常见错误：
  1. 先按 FDR<0.05 切一刀，再对"显著基因"做 Fisher 检验（**二分法丢掉了效应量的连续信息**，且对阈值敏感）；
  2. 忽略基因间的相关（共表达模块），把 100 个相关基因当成 100 个独立证据（**p 值严重膨胀**）。

---

## ② 背景原理

### 2.1 竞争性检验 vs 自含式检验

| 类型 | 零假设 | 例子 | 回答 |
|---|---|---|---|
| **自含式（self-contained）** | 该基因集内**没有任何基因**差异 | `limma::geneSetTest`（mean-rank） | "这个集合里有信号吗？" |
| **竞争性（competitive）** | 该基因集的差异程度**不高于集合外基因** | `limma::camera`、`fry`、GSEA | "这个集合比别的集合更差异吗？" |

原文用的是**竞争性**检验（camera 家族）。竞争性检验的零假设更贴近"通路特异性"，但更容易受"全局表达变化"影响。

### 2.2 camera 的统计量

设基因集 $S$ 有 $m$ 个基因，每个基因的 moderated t 为 $t_j$。

1. 集合内均值：
$$
\bar t=\frac{1}{m}\sum_{j\in S} t_j.
$$

2. 集合内 t 的样本方差：
$$
s_t^2=\frac{1}{m-1}\sum_{j\in S}(t_j-\bar t)^2.
$$

3. **基因间相关校正**：若集合内基因两两相关（平均相关 $\bar\rho$），则 $\bar t$ 的方差不是 $s_t^2/m$，而是
$$
\operatorname{Var}(\bar t)\approx\frac{s_t^2}{m}\bigl(1+(m-1)\bar\rho\bigr).
$$
定义方差膨胀因子 $\text{VIF}=1+(m-1)\bar\rho$。

4. 检验统计量：
$$
Z=\frac{\bar t}{\sqrt{s_t^2\,\text{VIF}/m}},
$$
在零假设下近似 $t$ 分布（自由度由设计矩阵决定；大样本近似正态）。

5. 双侧 p：
$$
p=2\,\Pr(T_{d}\ge|Z|).
$$

**直觉**：一个 200 个基因的集合，如果基因高度共表达（$\bar\rho=0.3$），VIF = 1+199×0.3 ≈ 60.7，有效样本量从 200 降到约 3.3——不校正就会得到荒谬的 p 值。

> **`zenith` 的实际实现是两样本形式**（与上式在集合远小于全体时近似等价）：
> $$
> \delta=\frac{G}{m_2}\left(\bar t_{\text{in}}-\bar t_{\text{all}}\right),\qquad
> \operatorname{Var}_{\text{pooled}}=\frac{(G-1)s_t^2-\delta^2 m m_2/G}{G-2},
> $$
> $$
> \mathrm{se}=\sqrt{\operatorname{Var}_{\text{pooled}}\left(\frac{\text{VIF}}{m}+\frac{1}{m_2}\right)},
> \qquad t=\frac{\delta}{\mathrm{se}},
> $$
> 其中 $m$ 为集合大小、$m_2=G-m$、$G$ 为全部基因数。输出里的 `delta`/`se` 就是这两个量。

### 2.3 camera 的两种实现

- `inter.gene.cor = 0.01`（limma 默认）：假设基因间相关很弱，直接用固定值。**保守**。
- `inter.gene.cor = NA`：从数据估计 $\bar\rho$（用基因集内基因的残差相关）。**更准，但更慢**。
- `use.ranks = TRUE`：基于秩的稳健版本（对离群基因不敏感）。

`zenith` 的默认是 `inter.gene.cor = 0.01`，并在此基础上适配 LMM（见 §2.4）。

### 2.4 zenith 对 LMM 的适配

`zenith::zenith()` 与 `limma::camera()` 的差别：

1. **输入不是表达矩阵，而是 dreamlet 的 t 统计量**（每基因一个 moderated t，已含权重与随机效应）。
2. **基因间相关从 LMM 残差估计**（而不是原始表达），避免随机效应造成的虚假相关。
3. **自由度**：LMM 下每基因的 $d_g$ 不同，zenith 用近似统一自由度或渐近正态。
4. 输出列（`zenith::zenith()` 的返回）：`NGenes`、`delta`（`= (G/m2)(集合内平均 t − 全体平均 t)`，见 §2.2 的两样本公式）、`se`、`p.less`、`p.greater`、`PValue`（双侧）、`FDR`、`Direction`；`dreamlet::zenith_gsa()` 再补 `assay`、`coef`、`Geneset`。**没有 `t` 列**；需要时自行派生 `t = delta/se`（作者代码正是这么做的）。

`dreamlet::zenith_gsa()` 是 `zenith::zenith_gsa()` 的 S4 方法，直接接受 `dreamletResult`：

```r
res_zenith <- zenith_gsa(fit, geneSets = go.gs, coefs = "DxAD", progressbar = FALSE)
```

### 2.5 基因集来源

| 来源 | 内容 | 本项目用途 | 版本 |
|---|---|---|---|
| **GO**（BP/CC/MF） | ~ 数万个条目，层级冗余 | Fig. 3G、6E 的"full spectrum"热图 | `GO.db` / `org.Hs.eg.db`（记录日期） |
| **MSigDB Hallmark** | 50 个精选集合，冗余低 | Fig. 3H 的 TNFα/NF-κB | MSigDB v2023.2+（记录） |
| **GO BP（p38/synapse）** | 层级条目；`GOBP_P38MAPK_CASCADE`、`GOBP_SYNAPSE_ASSEMBLY` | Fig. 6F 的 p38MAPK、Fig. 6E 的突触组装 | `GO.db` / `org.Hs.eg.db`（记录日期） |
| **MSigDB C2/C5** | 通路/GO | 可选扩展 | 同上 |

> ⚠️ **导师批注**：GO 的层级结构会导致同一信号在父子条目重复出现；原文 Fig. 3G 用 GO 做热图时只展示"每个细胞类型 top 5 + bottom 1"，就是为了避免冗余刷屏。MSigDB Hallmark 更适合做"可解释的通路结论"。

### 2.6 结果的方向与显著性

- `delta`（或 `t`）的**符号**：正 = 基因集整体上调；负 = 下调。
- `FDR`：BH 校正（跨基因集 × 细胞类型）。
- 原文 Fig. 3G/6E 的热图：行 = 基因集，列 = 细胞类型/疾病状态，颜色 = t 统计量，`*` = study-wide FDR<5%。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应操作 | 原文是否给参数 | 缺失项决策 |
|---|---|---|---|---|
| "Gene set analysis using the full spectrum of test statistics identified key immune response pathways activated in different cell types based on disease state." | 用全部 t 统计量做 GSA | `zenith_gsa()` | 未给基因集 | D9：GO + Hallmark 都跑 |
| "TNFα signaling via NF-κB showed activation in all COVID-19 disease states but not in patients with sepsis." | TNFα/NF-κB | Hallmark 集合 | 给了通路名 | — |
| "Within the TNFα signaling pathway genes showed different patterns based on disease state, with NFKBIA upregulation specific to COVID-19 patients, TNF activation only in mild and non-hospitalized COVID-19 patients, and CD83 upregulation specific to non-hospitalized COVID-19 patients." | 通路内基因异质性 | `plotGeneHeatmap` | 给了基因名 | — |
| "Gene set analysis … identifies up-regulation of synapse assembly and glutamate pathways in subsets of excitatory and inhibitory neurons." | 突触通路 | GO/Hallmark | 给了描述 | — |
| "microglia have a specific up-regulation of the p38MAPK cascade, which is involved in microglial inflammatory response." | p38MAPK | Hallmark/GO | 给了通路 | — |
| "This pathway includes DUSP10, which is one of the top upregulated genes in microglia." | DUSP10 | `topTable` + 热图 | 给了基因 | — |

---

## ④ 逐步可运行代码

### 4.1 加载基因集

```r
# code/07_00_genesets.R
suppressPackageStartupMessages({
  library(zenith); library(GSEABase); library(msigdbr); library(org.Hs.eg.db)
})

# (a) GO：to = "SYMBOL" 或 "ENSEMBL"；按需选 BP/CC/MF
go.gs <- get_GeneOntology("BP", to = "SYMBOL")
# 去掉过大的集合（作者代码：>1000 基因）
go.gs <- go.gs[sapply(go.gs, function(x) length(geneIds(x))) < 1000]

# (b) MSigDB Hallmark
# ⚠️ msigdbr >= 10 用 collection="H"；旧版用 category="H"（按当前包版本核对）
hall <- msigdbr(species = "Homo sapiens", collection = "H")
hall.gs <- split(hall$gene_symbol, hall$gs_name)
hall.gs <- lapply(hall.gs, unique)
hall.gsc <- GeneSetCollection(lapply(names(hall.gs), function(nm) {
  GeneSet(hall.gs[[nm]], setName = nm, geneIdType = SymbolIdentifier())
}))

saveRDS(list(go = go.gs, hallmark = hall.gsc), "processed/genesets.rds")

# 记录版本
writeLines(capture.output(sessionInfo()), "logs/07_sessionInfo.txt")
```

### 4.2 运行 zenith

```r
# code/07_01_zenith.R
# 用法：Rscript code/07_01_zenith.R psychad DxAD
suppressPackageStartupMessages({ library(dreamlet); library(zenith); library(GSEABase) })

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[1] else "psychad"
coef    <- if (length(args) >= 2) args[2] else "DxAD"

fit <- readRDS(sprintf("processed/fit_%s.rds", dataset))
gs  <- readRDS("processed/genesets.rds")

run_one <- function(gsc, tag) {
  res <- zenith_gsa(fit, geneSets = gsc, coefs = coef,
                    use.ranks = FALSE, n_genes_min = 10,
                    inter.gene.cor = 0.01, progressbar = FALSE)
  res$geneset_source <- tag
  write.table(res, sprintf("results/zenith_%s_%s_%s.tsv", dataset, coef, tag),
              sep = "\t", quote = FALSE, row.names = FALSE)
  res
}

res_go  <- run_one(gs$go, "GO_BP")
res_h   <- run_one(gs$hallmark, "Hallmark")
```

> 若 `zenith_gsa()` 报 "coefs are not found in dreamletResult"，先用 `coefNames(fit)` 确认；参数名是 **`coefs`**（不是 `coef`）。AD 的系数名是 `DxAD`，COVID 是 `SourceCOVID_MILD`，Kfoury 是 `Tumor_Involved` 等。

### 4.3 出图：原文风格热图

```r
# code/07_02_zenith_plots.R
suppressPackageStartupMessages({ library(zenith); library(dplyr); library(ggplot2) })

res <- read.delim("results/zenith_psychad_DxAD_Hallmark.tsv")

# (a) 每个细胞类型选 top 5 上调 + bottom 3 下调
p <- plotZenithResults(res, 5, 3)
ggsave("figures/Fig6E_zenith_AD.pdf", p, width = 14, height = 10)
ggsave("figures/Fig6E_zenith_AD.png", p, width = 14, height = 10, dpi = 300)

# (b) 指定通路热图（TNFα / p38MAPK / 突触）
# ⚠️ MSigDB Hallmark 只有 50 个集合，没有 P38_MAPK / SYNAPSE_ASSEMBLY；
# 这两个通路请用 GO BP 集合（见 §2.5）
sel <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
go_sel <- c("GOBP_P38MAPK_CASCADE", "GOBP_SYNAPSE_ASSEMBLY")
res_go <- read.delim("results/zenith_psychad_DxAD_GO_BP.tsv")
res_sel <- dplyr::bind_rows(
  res[res$Geneset %in% sel, ],
  res_go[res_go$Geneset %in% go_sel, ]
)
p2 <- ggplot(res_sel, aes(x = assay, y = Geneset, fill = delta)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "mean t")
ggsave("figures/07_selected_pathways_AD.png", p2, width = 9, height = 3, dpi = 300)
```

### 4.4 手动复算一个基因集的 camera 统计量（验证理解）

```r
# code/07_03_camera_manual.R
suppressPackageStartupMessages({ library(limma); library(GSEABase) })

fit  <- readRDS("processed/fit_psychad.rds")
gsc  <- readRDS("processed/genesets.rds")$hallmark
set_name <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
genes <- geneIds(gsc[[set_name]])

tt <- topTable(fit[["Micro_PVM"]], coef = "DxAD", number = Inf, sort.by = "none")
tt$ID <- rownames(tt)
sub <- tt[tt$ID %in% genes, ]

m <- nrow(sub)
tbar <- mean(sub$t)
s2   <- var(sub$t)
rho  <- 0.01                       # 与 zenith 默认一致
vif  <- 1 + (m - 1) * rho
Z    <- tbar / sqrt(s2 * vif / m)
df   <- min(fit[["Micro_PVM"]]$df.residual, na.rm = TRUE) + fit[["Micro_PVM"]]$df.prior
p    <- 2 * pt(-abs(Z), df = df)

cat(sprintf("set=%s m=%d tbar=%.3f s2=%.3f VIF=%.2f Z=%.3f p=%.3e\n",
            set_name, m, tbar, s2, vif, Z, p))

# 与 zenith 输出对照
z <- read.delim("results/zenith_psychad_DxAD_Hallmark.tsv")
print(z[z$Geneset == set_name & z$assay == "Micro_PVM", ])
```

> 数值不会完全相同：zenith 可能用不同的相关估计/自由度。**方向与数量级必须一致**；若差一个数量级，检查 `m`（基因匹配率）与 `inter.gene.cor`。

### 4.5 通路内基因热图（Fig. 3H 风格）

```r
# code/07_04_pathway_gene_heatmap.R
library(dreamlet); library(ggplot2); library(tidyr); library(dplyr)

fit <- readRDS("processed/fit_covid.rds")
# 从 MSigDB Hallmark TNFA 集合取基因（不要手写基因列表）
suppressPackageStartupMessages(library(msigdbr))
hall <- msigdbr(species = "Homo sapiens", collection = "H")
genes <- unique(hall$gene_symbol[hall$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"])
# 只保留在 cMono 结果中出现的基因
avail <- rownames(topTable(fit[["cMono"]], coef = "SourceCOVID_MILD", number = Inf, sort.by = "none"))
genes <- intersect(genes, avail)
stopifnot(length(genes) >= 10)

p <- plotGeneHeatmap(fit, coef = "SourceCOVID_MILD", genes = genes, assays = "cMono")
ggsave("figures/Fig3H_tnfa_genes.pdf", p, width = 10, height = 5)
ggsave("figures/Fig3H_tnfa_genes.png", p, width = 10, height = 5, dpi = 300)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| `geneSets` | GO BP + Hallmark | C2/C5、自定义 | 集合定义决定结论；必须在方法学写明 |
| `coefs` | 每个对比单独跑 | 多个 | 每个 coef 一行结果 |
| `use.ranks` | `FALSE` | `TRUE` | 秩检验更稳健但功效略低 |
| `inter.gene.cor` | `0.01`（默认） | `NA`（估计）、0.05 | 越大越保守；`NA` 最准但慢 |
| `n_genes_min` | 10 | 5–30 | 太小的集合不稳定 |
| `n_genes_max`（预处理） | 1000（GO） | — | 过大的集合信号被稀释 |
| FDR | BH，跨基因集×细胞类型 | 细胞类型内 | 影响"显著通路数" |
| 热图展示 | top 5 + bottom 3 | top N | 只影响可视化 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期（原文） | 允许波动 | Go/No-Go |
|---|---|---|---|
| COVID：TNFα/NF-κB 在 COVID 各状态激活、Sepsis 不激活 | 方向一致 | 定性 | ❌ |
| COVID：胆固醇稳态、IFN-γ、OXPHOS 在住院患者激活 | 方向一致 | 定性 | ⚠️ |
| Kfoury：MHC-II 与蛋白折叠 chaperone 在 tumor 上调 | 方向一致 | 定性 | ❌ |
| AD：microglia p38MAPK 上调（含 DUSP10） | 方向一致 | 定性 | ❌ |
| AD：神经元突触组装/谷氨酸通路上调 | 方向一致 | 定性 | ⚠️ |
| 显著基因集数 | 与原文同量级 | ±30% | ⚠️ |
| 基因匹配率（集合内基因在结果中的比例） | >70% | — | ❌ 低则查 ID 转换 |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| 所有通路都"显著" | 未校正基因间相关 / 用了自含式检验 | 用 camera/zenith 的竞争性检验；`inter.gene.cor=NA` |
| 通路方向与原文相反 | 系数符号/参考组反了 | 检查 `levels` 与 contrast |
| 基因匹配率极低 | 基因集是 SYMBOL，结果是 Ensembl | 统一 ID；打印交集大小 |
| `zenith_gsa` 报错 "geneSets must be GeneSetCollection" | 传了 list | 用 `GSEABase::GeneSetCollection()` 包裹 |
| 结果里 `NGenes` 很小 | `n_genes_min` 过滤 + 匹配失败 | 检查 ID 与集合版本 |
| GO 热图刷屏 | 层级冗余 | 只展示 top/bottom；或改用 Hallmark |
| 不同 MSigDB 版本结果不同 | 集合内容变化 | 记录版本与日期（L0 §5.5） |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 3G（COVID 通路热图） | `figures/Fig3G_gsa_covid.*` | 五个疾病状态 × 五个细胞亚型 |
| Fig. 3H（TNFα 通路基因热图） | `figures/Fig3H_tnfa_genes.*` | `plotGeneHeatmap()` |
| Fig. 4G（Kfoury 通路） | `figures/Fig4G_gsa_kfoury.*` | MHC-II / chaperone |
| Fig. 4H（chaperone 基因） | `figures/Fig4H_chaperone.*` | 基因级热图 |
| Fig. 6E（AD 通路） | `figures/Fig6E_zenith_AD.*` | `plotZenithResults()` |
| Fig. 6F（p38MAPK 基因） | `figures/Fig6F_p38.*` | DUSP10 等 |

> **交付物**：`results/zenith_<dataset>_<coef>_<source>.tsv`、`figures/Fig*G/H_*`、`logs/07_sessionInfo.txt`。完成后进入 `08_计算与统计性能基准复现.md`。
