# 08 代谢 CRISPR 筛选与靶向代谢组

> **目标文献**：França GS, Baron M, King BR, Bossowski JP, Bjornberg A, Pour M, Rao A, Patel AS, Misirlioglu S, Barkley D, Tang KH, Dolgalev I, Liberman DA, Avital G, Kuperwaser F, Chiodin M, Levine DA, Papagiannakopoulos T, Marusyk A, Lionnet T, Yanai I. *Cellular adaptation to cancer therapy along a resistance continuum.* **Nature (2024) 631:876–890.** DOI: 10.1038/s41586-024-07690-9
> **模块对应原文**：Fig. 5a–f；Extended Data Fig. 10a–h；Methods「Metabolic CRISPR–Cas9 genetic screen」「Data analysis of the metabolic CRISPR–Cas9 genetic screen」「Targeted metabolomics profiling」「Immunoblotting」
> **前置依赖**：`01`（GSE247686、ST003102）、`04`（模块 A–F）；本模块只产出"依赖清单与代谢组差异"，CB-839 联用表型由 09 号模块承接（不反向依赖）
> **复现难度 / 预计工时 / 计算资源**：★★★★☆（计算）/ 2–3 人日 / 32 GB RAM；MAGeCK 计算量中等

---

## ① 科学问题与逻辑位置

**问题**：

1. 沿连续谱，细胞是否逐步获得特定的代谢依赖？
2. 这些依赖是否对应模块 E/F 的表达程序（NRF2/PPP/核苷酸/铁稳态）？
3. 代谢组是否支持核苷酸池升高？
4. 依赖是否可被药物组合利用（奥拉帕利 + CB-839）？

**逻辑位置**：04 号模块说"表达了什么"，本模块说"细胞**离不开**什么"。CRISPR 负筛提供功能遗传学证据，代谢组提供代谢物证据，09 号模块提供药理学验证。

**上游**：`raw/geo/GSE247686/*`、`raw/metabolomics/ST003102/*`、`results/04_gene_modules.csv`。
**下游**：`results/08_crispr_gene_scores.csv`、`results/08_metabolomics_diff.csv`、`figures/Fig5b_crispr_heatmap.*`。

---

## ② 背景原理

### 2.1 负选择 CRISPR 筛选

- 文库：每个基因 ~10 条 sgRNA；慢病毒低 MOI（~0.3）保证大多数细胞只整合一个 sgRNA。
- 嘌呤霉素选择后取 **t0**（初始）与 **end**（~10 个群体倍增后）。
- 若某基因被敲除后细胞增殖受损，其 sgRNA 在 end 中**减少**（depletion）。
- **MAGeCK-MLE** 用最大似然估计同时考虑所有 sgRNA，输出基因-level β score（负值 = 敲除有害）。
- 比较不同样本的 β score 差异，可识别"条件特异依赖"。

### 2.2 为什么用负选择而不是正选择

代谢依赖通常是**必需性**（loss-of-function 降低适应度），因此负选择更直接。原文还比较了"有药/无药"（CD、T320D），以区分"基础必需"与"药物暴露下必需"。

### 2.3 高维 CRISPR 的噪声

- sgRNA 效率差异；
- 文库覆盖不均；
- 细胞周期/倍性；
- 脱靶；
- 群体瓶颈。

因此原文要求：负 β、P<0.05、FDR<0.25；跨样本比较时还要求 |score difference| > 1 SD。**不要**只按 P 值挑基因。

### 2.4 靶向代谢组

- ZIC-pHILIC 色谱分离极性代谢物；Q Exactive HF 高分辨质谱；正负模式切换。
- 峰提取基于保留时间 + 精确质量（±5 ppm）；S/N ≥3、下限 10,000。
- 组间比较用 t 检验（原文 SciPy 1.1.0），P<0.05。
- 关键陷阱：代谢物名称/离子加合物、批次、空白校正、缺失值填补（用阈值）。

### 2.5 CB-839 的机制假设

奥拉帕利 → DNA 损伤/ROS → NRF2 抗氧化程序（GSH、PPP）→ 核苷酸合成前体；GLS 抑制谷氨酰胺→谷氨酸，削弱 GSH 与 TCA/核苷酸前体。两条通路同时被压制 → 协同。Bliss 独立模型用于量化协同：

$$E_{AB} = E_A + E_B - E_A E_B,\quad S = E_{obs} - E_{AB}$$

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 操作 | 原文参数 | 复现决策 |
|---|---|---|---|---|
| "metabolic gene-knockout library (2,978 genes, 29,790 gRNAs)" | 代谢文库 | 文库 | 2,978/29,790 | 用公开序列 |
| "Kuramochi C, T10 and T320 (CD44hiCD24low, sorted)" | 分选后感染 | 细胞 | — | 方案 |
| "multiplicity of infection of ~0.3" | MOI 0.3 | 感染 | 0.3 | 固定 |
| "puromycin (4 μg/ml) 48 h after infection" | 选择 | — | 4 μg/ml | 固定 |
| "cultured for 10 population doublings" | 10 倍增 | 筛选窗口 | 10 | 固定 |
| "C in 1 μM and T320 in 40 μM" | 加药条件 | CD/T320D | 1/40 μM | 固定 |
| "MAGeCK-count ... MAGeCK-MLE ... MAGeCK 0.5.8" | 计数与打分 | MAGeCK | 0.5.8 | 固定 |
| "negative β scores, P < 0.05 and FDR < 0.25" | 依赖判据 | — | — | 固定 |
| "negative score differences higher than one standard deviation" | 样本间差异 | — | 1 SD | 固定 |
| "T320D included comparisons with T320 and CD" | 双重比较 | — | — | 固定 |
| "fgsea ... P < 0.05 and FDR < 0.25" | 通路 | fgsea | — | 固定 |
| "ZIC-pHILIC ... Q Exactive HF ... 5 ppm ... S/N 3× ... floor 10,000" | 代谢组 | LC-MS | — | 固定 |
| "SciPy 1.1.0 ... t-tests ... P < 0.05" | 统计 | t 检验 | — | 固定 |
| "CB-839 (50 nM)" | GLS 抑制剂 | 联用 | 50 nM | 固定 |
| "sixfold decrease in IC50" | 最适应系 | 联合 | — | 09 号模块 |
| "Bliss synergy score" | 协同 | — | — | 固定 |

**原文未报告**：文库的 sgRNA 序列版本/Addgene 编号；测序深度（20–40 M paired reads 已给）；每个样本的细胞数；代谢组的完整归一化流程；缺失值处理细节。均登记 D 表。

---

## ④ 逐步可运行代码

### 4.1 整理 MAGeCK 输入

```r
# code/08_01_prepare_counts.R
suppressPackageStartupMessages({library(data.table); library(dplyr)})
cnt <- fread("raw/geo/GSE247686/GSE247686_kuramochi_C_T10_T320_sgrna_counts.tsv.gz",
             data.table = FALSE)
lib <- fread("raw/geo/GSE247686/GSE247686_Large_metabolic_library_sgRNA_sequence_human.tsv.gz",
             data.table = FALSE)

cat("counts columns:\n"); print(colnames(cnt))
cat("library columns:\n"); print(colnames(lib))

# 实测 counts 列：sgRNA, Gene, C_t0, C_end, C_drug_end, T10_t0, T10_end,
#                  T40_t0, T40_end, T320_t0, T320_end, T320_drug_end
sgrna_col <- "sgRNA"; gene_col <- "Gene"
stopifnot(all(c(sgrna_col, gene_col) %in% colnames(cnt)))
sample_cols <- setdiff(colnames(cnt), c(sgrna_col, gene_col))
print(sample_cols)

fwrite(cnt[, c(sgrna_col, gene_col, sample_cols)],
       "processed/08_mageck_counts.tsv", sep = "\t")

# MAGeCK-MLE 设计矩阵：第一列必须叫 Samples；baseline 全为 1；
# 每个条件列 1 = 该样本属于该条件；t0 样本所有条件列为 0。
design <- data.frame(
  Samples = sample_cols,
  baseline = 1,
  C     = as.integer(sample_cols %in% "C_end"),
  CD    = as.integer(sample_cols %in% "C_drug_end"),
  T10   = as.integer(sample_cols %in% "T10_end"),
  T40   = as.integer(sample_cols %in% "T40_end"),      # 非主分析，敏感性用（D10）
  T320  = as.integer(sample_cols %in% "T320_end"),
  T320D = as.integer(sample_cols %in% "T320_drug_end"),
  stringsAsFactors = FALSE
)
write.table(design, "processed/08_mageck_design.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(design)
```

> **实测说明**：GEO 子系列登记 8 个样本，但 counts 表含 **10 个样本列**（多出 T40_t0/T40_end）。主分析按原文只用 C/CD/T10/T320/T320D；T40 仅作敏感性（D10）。`baseline` 列必须全为 1，否则 MAGeCK 无法估计共同基线。

### 4.2 运行 MAGeCK-MLE

```bash
#!/usr/bin/env bash
# code/08_02_mageck_mle.sh
set -euo pipefail
OUT=results/08_mageck
mkdir -p "$OUT"

mageck mle \
  -k processed/08_mageck_counts.tsv \
  -d processed/08_mageck_design.txt \
  -n "$OUT/mageck_mle" \
  --norm-method median \
  --threads 8

# 输出：mageck_mle.gene_summary.txt / .sgrna_summary.txt
head -2 "$OUT/mageck_mle.gene_summary.txt"
```

### 4.3 基因打分与样本间差异

```r
# code/08_03_gene_scores.R
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})
gs <- read_tsv("results/08_mageck/mageck_mle.gene_summary.txt", show_col_types = FALSE)
# 列名形如：Gene, C|beta, C|z, C|p-value, C|fdr, ...
beta_cols <- grep("\\|beta$", colnames(gs), value = TRUE)
fdr_cols  <- grep("\\|fdr$", colnames(gs), value = TRUE)
p_cols    <- grep("\\|p-value$", colnames(gs), value = TRUE)

long <- gs %>%
  select(Gene, all_of(beta_cols), all_of(p_cols), all_of(fdr_cols)) %>%
  pivot_longer(-Gene, names_to = c("sample", ".value"), names_pattern = "(.*)\\|(.*)") %>%
  rename(beta = beta, p = `p-value`, fdr = fdr)

write.csv(long, "results/08_crispr_gene_scores.csv", row.names = FALSE)

# 样本 vs 参考的 score difference
ref <- "C"
comp <- long %>%
  filter(sample %in% c("C", "CD", "T10", "T320", "T320D")) %>%
  select(Gene, sample, beta, p, fdr) %>%
  pivot_wider(names_from = sample, values_from = c(beta, p, fdr))

for (s in c("CD", "T10", "T320", "T320D")) {
  comp[[paste0("delta_", s)]] <- comp[[paste0("beta_", s)]] - comp[[paste0("beta_", ref)]]
}

# 依赖判据：负 beta、P<0.05、FDR<0.25、|delta| > 1 SD
sd_delta <- sd(comp$delta_T320, na.rm = TRUE)
dep <- comp %>%
  filter(beta_T320 < 0, p_T320 < 0.05, fdr_T320 < 0.25,
         delta_T320 < -sd_delta)
write.csv(dep, "results/08_crispr_depleted_T320.csv", row.names = FALSE)
cat("depleted genes (T320 vs C):", nrow(dep), "\n")
```

### 4.4 通路富集

```r
# code/08_04_crispr_gsea.R
suppressPackageStartupMessages({
  library(fgsea); library(msigdbr); library(dplyr); library(tidyr)
})
# 从 4.3 的 long 表重建条件宽表（脚本独立运行，不依赖内存中的 comp）
long <- read.csv("results/08_crispr_gene_scores.csv")
comp <- long %>%
  select(Gene, sample, beta, p, fdr) %>%
  pivot_wider(names_from = sample, values_from = c(beta, p, fdr))
comp$delta_T320 <- comp$beta_T320 - comp$beta_C

ranked <- comp %>%
  filter(!is.na(delta_T320)) %>%
  mutate(rank = -delta_T320) %>%       # 越负越靠前
  arrange(desc(rank))
stats <- setNames(ranked$rank, ranked$Gene)

reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
pathways <- bind_rows(reactome, kegg) %>% split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

set.seed(20261018)
res <- fgsea(pathways, stats, minSize = 10, maxSize = 500, nperm = 10000)
res <- res %>% filter(pval < 0.05, padj < 0.25) %>% arrange(pval)
write.csv(res, "results/08_crispr_gsea.csv", row.names = FALSE)
print(res[, c("pathway","NES","pval","padj")])
```

### 4.5 代谢组处理

```r
# code/08_05_metabolites.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ComplexHeatmap)
})
# 输入格式：行=代谢物，列=样本（C_1..C_3, T320_1..T320_3, C_1D.., T320_1D..）
met <- read.delim("raw/metabolomics/ST003102/metabolite_intensity.tsv",
                  row.names = 1, check.names = FALSE)

# log2 + 中位数归一化
logm <- log2(met + 1)
norm <- sweep(logm, 2, apply(logm, 2, median, na.rm = TRUE), "-")

# 分组
grp <- ifelse(grepl("^C", colnames(norm)), "C", "T320")
dose <- ifelse(grepl("D$", colnames(norm)), "drug", "nodrug")

res <- lapply(unique(dose), function(d) {
  idx <- dose == d
  x <- norm[, idx, drop = FALSE]
  g <- grp[idx]
  p <- apply(x, 1, function(v) t.test(v[g == "C"], v[g == "T320"])$p.value)
  fc <- apply(x, 1, function(v) mean(v[g == "T320"]) - mean(v[g == "C"]))
  data.frame(metabolite = rownames(x), comparison = d,
             log2FC = fc, p = p)
}) %>% bind_rows() %>% mutate(p_adj = p.adjust(p, method = "BH"))

write.csv(res, "results/08_metabolomics_diff.csv", row.names = FALSE)
sig <- res %>% filter(p < 0.05) %>% arrange(p)
print(head(sig, 20))
```

### 4.6 CB-839 联用曲线（Bliss 计算在 09 号模块）

> 本模块只提供"为什么联用"的代谢依据；**IC50 与 Bliss 的正式计算在 09 号模块**（`code/09_05_combination.R`），避免两处实现不一致。下图仅用于展示联用曲线。

```r
# code/08_06_combination_plot.R
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
# 输入：results/09_viability_combination.csv
# 列：line, olaparib_uM, cb839_nM, replicate, viability
d <- read.csv("results/09_viability_combination.csv")
avg <- d %>% group_by(line, olaparib_uM, cb839_nM) %>%
  summarise(viability = mean(viability), .groups = "drop")

p <- ggplot(avg, aes(olaparib_uM, viability, colour = factor(cb839_nM))) +
  geom_point() + geom_line() + scale_x_log10() +
  labs(x = "Olaparib (μM)", y = "Viability (%)", colour = "CB-839 (nM)") +
  theme_classic()
ggsave("figures/Fig5d_cb839_curves.pdf", p, width = 5, height = 4)
ggsave("figures/Fig5d_cb839_curves.png", p, width = 5, height = 4, dpi = 300)
```

---

## ⑤ 参数解释与调参空间

| 参数 | 原文 | 本项目 | 改变后果 |
|---|---|---|---|
| MOI | 0.3 | 0.3 | 过高导致多 sgRNA/细胞 |
| 嘌呤霉素 | 4 μg/ml | 4 | 过低选择不纯 |
| 倍增数 | 10 | 10 | 过少信号弱，过多饱和 |
| 文库覆盖 | ~10 sgRNA/基因 | 同 | 覆盖不足增加噪声 |
| 归一化 | median | median | 其他方法改变 β |
| 依赖阈值 | β<0、P<0.05、FDR<0.25 | 同 | 放宽增加假阳性 |
| 样本差异 | >1 SD | 同 | 过严漏掉弱依赖 |
| GSEA | P<0.05、FDR<0.25 | 同 | 同上 |
| 代谢组 t 检验 | P<0.05 | 同 + BH | 原文未说明多重校正；双口径 |
| CB-839 | 50 nM | 50 nM | 浓度过高会产生单药毒性 |
| Bliss | 独立模型 | 同 | 只适用于独立作用假设 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 原文/预期 | 判据 |
|---|---|---|
| 文库基因数 | 2,978 | 文件核对 |
| sgRNA 数 | 29,790 | 文件核对 |
| 样本数 | GEO 登记 8 个；counts 表 10 个样本列（含 T40） | 列名核对；主分析排除 T40 |
| mapping rate | 75–80% | `mageck count` 日志 |
| 依赖基因数 | T10 < T320 | 单调趋势 |
| 共同依赖 | TXNRD1、PRDX1、PGD、PGK1、ALDOA、ENO1、STEAP3、IMPDH2、GUK1、DTYMK | ≥70% 命中 |
| T320/T320D 额外 | TYMS、CAD、CTPS1、FDX1、FXN | 方向一致 |
| T320D 特异 | COX5B、CYCS、ACAD9 | 方向一致 |
| 代谢组 | 核苷酸/前体升高 | ≥3 个 P<0.05 |
| CB-839 | IC50 降约 6 倍 | ≥3 倍，方向一致 |
| 单药 CB-839 | 轻微敏感 | 与联合对比 |
| 重复性 | 两次分析 | 关键基因排名稳定 |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| MAGeCK 报设计矩阵共线 | C/CD、T320/T320D 同时为 1 | 按互斥重编码；`baseline` 必须全为 1（4.1） |
| MAGeCK 报列名不匹配 | 设计矩阵第一列不是 `Samples` | 改为 `Samples`；样本名与 counts 完全一致 |
| mapping rate 低 | 文库版本不匹配 | 用公开 library 序列；检查 sgRNA ID |
| 所有基因都"必需" | 归一化/选择压力异常 | 检查 t0 与 end 的文库复杂度 |
| 关键基因未检出 | sgRNA 覆盖不足 | 报告覆盖；不要强行解释 |
| 代谢物缺失值多 | 低于检测限 | 用阈值填补；敏感性分析 |
| 代谢组批次效应 | 样本处理顺序 | 中位数归一化；QC 样本 |
| CB-839 单药毒性大 | 浓度过高 | 滴定 10–100 nM |
| Bliss 分数虚高 | 单药曲线未覆盖 | 全剂量矩阵 |
| 免疫印迹定量不准 | 曝光饱和 | 线性范围；ImageJ 校正 |
| 结果不可重复 | 细胞代次/分选 | 固定分选门；记录代次 |

---

## ⑧ 与文献图表的对应关系

| 原文图 | 内容 | 本项目产物 | 说明 |
|---|---|---|---|
| Fig. 5a | 筛选设计 | `figures/Fig5a_screen_scheme.pdf/png` | 方案图 |
| Fig. 5b | 差异依赖热图 | `figures/Fig5b_crispr_heatmap.pdf/png` | 264 基因（原文） |
| Fig. 5c | 通路 GSEA | `results/08_crispr_gsea.csv` | 核苷酸/PPP/OXPHOS/铁 |
| Fig. 5d | CB-839 机制示意 | `figures/Fig5d_cb839_scheme.pdf/png` | 方案图 |
| Fig. 5e | IC50 ratio | `figures/Fig5e_cb839_ic50.pdf/png` | 6 倍 |
| Extended Data Fig. 10a | 排名图 | `figures/FigED10a_rank.pdf/png` | CD/T10/T320/T320D |
| Extended Data Fig. 10b | 重叠韦恩/超几何 | `results/08_depleted_overlap.csv` | |
| Extended Data Fig. 10c | 代谢组热图 | `figures/FigED10c_metabolites.pdf/png` | |
| Extended Data Fig. 10d–f | 联合曲线/IC50/单药 | `figures/FigED10d_f_*.pdf/png` | 09 号模块 |
| Extended Data Fig. 10g | Bliss 热图 | `figures/FigED10g_bliss.pdf/png` | |
| Extended Data Fig. 10h | NRF2/HO-1 WB | `figures/FigED10h_wb.pdf/png` | 09 号模块 |

> **交付判定**：MAGeCK 运行成功、关键依赖基因方向一致、代谢组核苷酸升高、CB-839 联用 IC50 下降 ≥3 倍，方可进入 09 号模块。
