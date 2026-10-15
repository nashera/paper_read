# 03 细胞比例 Ro/e 富集分析与代谢通路 GSEA

> **前置依赖**：02。**对应原文**：Fig. 2F–H、Fig. S1C–D；Suppl. Table 3。
> **产出**：三组与 eMASH/aMASH 的 Ro/e 矩阵与统计检验、各细胞簇代谢通路 GSEA 热图。

---

## 1. Ro/e 的统计定义与原理

**Ro/e（observed/expected ratio）** 用于回答"某细胞类型在某疾病组中是富集还是耗竭"，同时消除各样本细胞总数不同的偏差。

设细胞类型 $c$、疾病组 $g$、样本 $s$：

- 观察值 $O_{c,s}$ = 样本 $s$ 中类型 $c$ 的细胞数占该样本总细胞数的比例（或计数）；
- 期望值 $E_{c,s} = n_s \times p_c$，其中 $n_s$ 为样本 $s$ 的细胞总数，$p_c$ 为类型 $c$ 在全队列的占比——即假设"细胞类型组成与疾病无关"时按边际分布计算的期望。

常用实现（原文 Fig. 2E–G 的做法）：

$$\mathrm{Ro/e}_{c,g} = \frac{\sum_{s \in g} O_{c,s}}{\sum_{s \in g} E_{c,s}}$$

- $\mathrm{Ro/e} > 1$：该组中该类型富集；$< 1$：耗竭。
- **统计检验**：对每个 $(c, g)$，在样本水平用 Wilcoxon rank-sum 检验比较 Ro/e（或比例）与其余组；或用卡方/置换检验在细胞水平检验"类型×组"的独立性。原文统计约定为 Wilcoxon rank-sum（Methods: Statistical analysis）。

> 直觉：Ro/e 本质是"条件概率之比"——$P(c \mid g) / P(c)$。它把"某组样本测序深度大"这种技术因素除掉，只留组成差异。

---

## 2. 代码（R）

```r
# code/03_roe_gsea.R
library(Seurat); library(ggplot2); library(dplyr); library(reshape2)
so <- readRDS("processed/integrated_harmony.rds")

prop <- table(so$celltype, so$sample) / table(so$sample)      # 样本×类型 比例矩阵
prop <- as.data.frame.matrix(prop)
meta <- so@meta.data |> distinct(sample, disease, dataset, mash_sub) |> as.data.frame()

# --- Ro/e：以全队列类型占比为期望 ---
p_all <- colMeans(prop)                                        # 类型全局占比（按样本平均）
roe <- prop / rep(p_all, each = nrow(prop))
roe$disease <- meta$disease[match(rownames(roe), meta$sample)]

roe_long <- melt(roe, id.vars = "disease", variable.name = "celltype")
roe_sum <- roe_long |> group_by(disease, celltype) |>
  summarise(RoE = mean(value), sd = sd(value), n = n(), .groups = "drop")

# --- 组间检验（Wilcoxon rank-sum，原文统计约定）---
pw <- list()
for (ct in colnames(prop)) {
  for (pr in list(c("MASH","MASLD"), c("MASH","Healthy"), c("MASLD","Healthy"))) {
    a <- roe[roe$disease == pr[1], ct]; b <- roe[roe$disease == pr[2], ct]
    p <- wilcox.test(a, b)$p.value
    pw[[length(pw)+1]] <- data.frame(celltype = ct, group1 = pr[1], group2 = pr[2],
                                     p = p, RoE1 = mean(a), RoE2 = mean(b))
  }
}
pw <- bind_rows(pw); pw$p_adj <- p.adjust(pw$p, "BH")
fwrite(pw, "results/Fig2F_roe_tests.csv")

# --- 绘图：Fig.2F（三组）---
ggplot(roe_sum, aes(celltype, disease, fill = RoE)) + geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "blue", mid = "white", high = "red") +
  theme_minimal() + coord_flip()
ggsave("figures/Fig2F_RoE_heatmap.pdf", width = 8, height = 5)

# --- Fig.2G：eMASH vs aMASH（仅 GSE289173 有分期）---
sub <- meta[!is.na(meta$mash_sub), ]
roe2 <- prop[sub$sample, , drop = FALSE] / rep(p_all, each = nrow(sub))
roe2$mash_sub <- sub$mash_sub
# 同上汇总 + wilcox，重点看 HSC 在 aMASH > eMASH
```

**核对点**：mono/mac 与 HSC 在 MASH 的 Ro/e 显著 > MASLD/Healthy；aHSC（04 文档细分后）在 aMASH > eMASH。

---

## 3. 各细胞簇代谢通路 GSEA（Fig. 2H）

### 3.1 GSEA 原理（预排序版, pre-ranked）

对每个细胞簇 $k$：
1. 用 `FindMarkers`（wilcox）得到簇 $k$ vs 其余所有细胞的 $\log_2\mathrm{FC}$ 与 $p$ 值；
2. 排序统计量 $r_g = \mathrm{sign}(\log_2\mathrm{FC}) \times (-\log_{10} p)$（或直接用 $\log_2\mathrm{FC}$）；
3. 对基因集 $S$（代谢通路，如 KEGG/Reactome 或原文的代谢基因集），沿排序基因列表 $L$ 从第 1 位走到第 $j$ 位，累积命中得分与未命中罚分：

$$\mathrm{ES}(S) = \max_{1 \le j \le N}\left| \sum_{g \in S,\, g \in L_{\le j}} \frac{|r_g|^p}{N_R} - \sum_{g \notin S,\, g \in L_{\le j}} \frac{1}{N - N_H} \right|$$

其中 $N_R = \sum_{g \in S}|r_g|^p$（$p=1$ 时为加权秩和相关），$N_H$ 为 $S$ 内基因数，$N$ 为总基因数；ES 取行走曲线偏离 0 的最大绝对值（GSEA 报告带符号版本）；
4. 以基因集大小匹配的置换（或基因置换）得到归一化富集分数 NES 与 FDR $q$ 值。

**单细胞场景的替代**：AUCell/UGSEA 亦可，但原文 Fig. 2H 为"簇 × 通路"热图，用 cluster-level pre-ranked GSEA 最直接。

### 3.2 代码

```r
library(clusterProfiler); library(org.Hs.eg.db)
# 代谢通路基因集：MSigDB Hallmark 之外，用 KEGG/Reactome 的代谢子集，
# 或自定义代谢基因集（原文 Fig.2H 为 metabolism-related pathways）
gmt <- read.gmt("data/gmt/kegg_metabolism.gmt")   # 或 c2.cp.reactome 过滤代谢词

res <- list()
for (ct in unique(so$celltype)) {
  mk <- FindMarkers(so, ident.1 = ct, group.by = "celltype",
                    logfc.threshold = 0.1, min.pct = 0.1, test.use = "wilcox")
  mk$gene <- rownames(mk)
  rank <- sort(mk$avg_log2FC, decreasing = TRUE)
  g <- GSEA(rank, TERM2GENE = gmt, pvalueCutoff = 0.05, seed = TRUE)
  res[[ct]] <- data.frame(celltype = ct, pathway = g$Description,
                          NES = g$NES, q = g$p.adjust)
}
res <- bind_rows(res)
# 热图：NES 矩阵（行=通路 top union，列=细胞类型），显著者标星
mat <- acast(res, pathway ~ celltype, value.var = "NES")
pheatmap(mat, na.col = "grey90", color = colorRampPalette(c("blue","white","red"))(100))
ggsave("figures/Fig2H_metabolism_GSEA.pdf", width = 10, height = 9)
```

**核对点**：肝细胞簇应出现显著代谢重编程（原文：hepatocytes 代谢通路显著变化），即肝细胞列的 |NES| 分布明显宽于其他列。

---

## 4. 延伸：把 Ro/e 推广到亚群层级

04/05 文档的亚群 Ro/e（HSC 5 亚型、髓系 7 亚群）复用本章同一函数，仅把 `celltype` 换成亚群标签、把分组换成 (MASLD, MASH) 或 (eMASH, aMASH)。**保持同一实现**可保证全文 Ro/e 口径一致（审稿人常查这一点）。

## 5. 陷阱

1. **样本量不均衡**：MASH 66 vs MASLD 14 vs Healthy 35，wilcox 功效不对称；报告效应量（Ro/e 差值）而不只看 p。
2. **比例和为 1 的组成数据（compositional）问题**：某类型升高会机械性压低其他类型比例。Ro/e 只描述相对组成；若要"绝对浸润量"结论需空间/流式等正交证据（原文用空间共定位补强）。
3. **多重检验**：类型×组×配对数较多，务必 BH 校正（代码已含）。
