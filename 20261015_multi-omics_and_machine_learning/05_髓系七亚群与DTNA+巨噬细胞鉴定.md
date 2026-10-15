# 05 髓系七亚群细分、DTNA⁺ 巨噬细胞鉴定与功能特征

> **前置依赖**：02（03 的 Ro/e 函数复用）。**对应原文**：Fig. 4 全图、Fig. 5A–D；Suppl. Table 6。
> **产出**：髓系 7 亚群注释、DTNA⁺ Mac 的 MASH 富集证据、代谢/炎症通路 z-score 与 Hallmark/irGSEA、髓系 CytoTRACE 与 State 轨迹。

---

## 1. 生物学背景与亚群定义

肝髓系 compartment 含组织驻留 Kupffer 细胞（KC）与单核来源巨噬细胞/DC。原文分辨 7 群：

| 亚群 | marker（原文 Fig.4B） | 身份 |
|---|---|---|
| VCAM1⁺ KC | VCAM1, CD163, MARCO | 驻留 KC 亚型 1 |
| SAA1⁺ KC | SAA1 | 驻留 KC 亚型 2（急性期反应样） |
| **DTNA⁺ Mac** | **DTNA, ITGAX** | 本研究核心新亚群 |
| VCAN⁺ Mon | VCAN, FCN1 | 经典单核（新浸润） |
| mDC | （髓系 DC marker，如 CLEC9A-/CD1C+ 谱） | 髓系树突细胞 |
| pDC | （如 LILRA4, IL3RA） | 浆细胞样 DC |
| MKI67⁺ Mac | MKI67, TOP2A | 增殖期巨噬 |

**DTNA**（dystrobrevin alpha）并非常规免疫 marker，这是数据驱动发现的亚群标志；ITGAX=CD11c，与既往"CD11c⁺ 巨噬驱动 MASH 纤维化"文献衔接（原文引 Itoh 2017）。

---

## 2. 子聚类与注释

```r
# code/05_myeloid.R
so <- readRDS("processed/integrated_harmony.rds")
mye <- subset(so, celltype == "MonoMac")        # n ≈ 37,907
mye <- NormalizeData(mye) |> FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |> RunPCA(npcs = 30) |>
       RunHarmony(group.by.vars = "dataset", dims.use = 1:20) |>
       FindNeighbors(reduction = "harmony", dims = 1:20) |>
       FindClusters(resolution = 0.6) |>
       RunUMAP(reduction = "harmony", dims = 1:20, seed.use = 20260101)
DimPlot(mye, label = TRUE)                       # Fig.4A：期望 7 群

mye_mk <- list(
  `VCAM1+ KC`  = c("VCAM1","CD163","MARCO","TIMD4"),
  `SAA1+ KC`   = c("SAA1","CD163"),
  `DTNA+ Mac`  = c("DTNA","ITGAX"),
  `VCAN+ Mon`  = c("VCAN","FCN1","S100A8"),
  mDC          = c("CD1C","CLEC10A"),
  pDC          = c("LILRA4","IL3RA","CLEC4C"),
  `MKI67+ Mac` = c("MKI67","TOP2A"))
# 02 文档打分-最大法注释；DotPlot 复核（Fig.4B）
DotPlot(mye, features = unique(unlist(mye_mk))) + RotatedAxis()
ggsave("figures/Fig4B_myeloid_dotplot.pdf")
saveRDS(mye, "processed/myeloid_annotated.rds")
```

---

## 3. Ro/e：DTNA⁺ Mac 在 MASH 富集（Fig. 4C）

复用 03 的 Ro/e 实现，分组 = (MASLD, MASH)。期望：DTNA⁺ Mac 的 Ro/e(MASH) 显著升高，VCAM1⁺ KC 下降。

---

## 4. 代谢/炎症通路 z-score 热图与箱线图（Fig. 4D–E）

### 4.1 原理

对每个预定义通路基因集 $S$（glycolysis、hypoxia、lipid mediator、M1/M2 polarization、pro-inflammatory、TCA cycle、interferon 等），用 `AddModuleScore`（Seurat）计算单细胞得分：对 $S$ 内基因取平均表达，减去"表达水平匹配的背景基因 bin"的平均（控制基因表达量对得分的影响），得到 per-cell score：

$$\mathrm{score}_i(S) = \frac{1}{|S|}\sum_{g \in S} x_{gi} \;-\; \frac{1}{|B|}\sum_{g \in B} x_{gi}$$

其中背景集 $B$ 由与 $S$ 内基因平均表达水平相近的基因组成（Seurat 将全部基因按平均表达分成 ~24 个 bin，对 $S$ 中每个基因从同 bin 随机抽 100 个组成 $B$）。随后：

- **热图（Fig.4D）**：对 (亚群 × 通路) 的平均 score 做 **z-score 标准化** $z = (s - \bar{s})/\mathrm{sd}(s)$（按行或按列，原文为 z-scored）；
- **箱线图（Fig.4E）**：仅在 DTNA⁺ Mac 内比较 MASLD vs MASH 的各通路 score，Wilcoxon 检验（原文 $p<0.05 \sim p<0.001$）。

通路基因集来源建议：MSigDB Hallmark（HALLMARK_HYPOXIA、HALLMARK_GLYCOLYSIS、HALLMARK_TNFA_SIGNALING_VIA_NFKB、HALLMARK_IL6_JAK_STAT3、HALLMARK_INTERFERON_GAMMA_RESPONSE、HALLMARK_COAGULATION、HALLMARK_CHOLESTEROL_HOMEOSTASIS）+ M1/M2 极化自定义集（M2: MRC1, CD163, ARG1, IL10…; M1: IL1B, TNF, CXCL9/10/11…）。

```r
hall <- read.gmt("data/gmt/h.all.v2023.1.Hs.symbols.gmt")
sets <- hall[hall$term %in% c("HALLMARK_HYPOXIA","HALLMARK_GLYCOLYSIS",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_IL6_JAK_STAT3",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_COAGULATION",
        "HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_ANGIOGENESIS"), ]
m2set <- list(M2 = c("MRC1","CD163","ARG1","IL10","TGFB1","CCL24"))
m1set <- list(M1 = c("IL1B","TNF","CXCL9","CXCL10","CXCL11","NOS2"))

mye <- AddModuleScore(mye, features = split(sets$gene, sets$term), name = "hs_")
mye <- AddModuleScore(mye, features = c(m1set, m2set), name = "pol_")

# Fig.4D 热图：亚群×通路 平均 score 的 z-score
mat <- as.data.frame.matrix(tapply(...) )   # 行=通路, 列=亚群, 值=mean score
pmat <- t(scale(t(mat)))                     # 行方向 z-score
pheatmap(pmat, color = colorRampPalette(c("blue","white","red"))(100))
ggsave("figures/Fig4D_pathway_zscore.pdf")

# Fig.4E：DTNA+ Mac 内 MASH vs MASLD 箱线 + wilcox
sub <- subset(mye, subtype == "DTNA+ Mac")
for (s in score_cols) {
  p <- wilcox.test(sub[[s]][sub$disease=="MASH"], sub[[s]][sub$disease=="MASLD"])$p.value
}
ggplot(sub, aes(disease, .data[[s]], fill = disease)) + geom_boxplot() +
  stat_compare_means(method = "wilcox.test")
```

**核对**：MASH 的 DTNA⁺ Mac 中 glycolysis、hypoxia、M2、pro-inflammatory、TCA、interferon 均显著升高（原文 Fig.4E 方向）。

---

## 5. Hallmark GSEA 热图（Fig. 4F）与 irGSEA（Fig. 4G）

```r
# Fig.4F：7 个髓系亚群各自的 Hallmark 富集热图（NES 矩阵）
# 做法同 03 的 cluster-level pre-ranked GSEA，群别换成髓系亚群
# 期望：DTNA+ Mac 列 TNFA-NFKB、IL6-JAK-STAT3、COAGULATION、CHOLESTEROL 高

# Fig.4G：irGSEA（单细胞水平 GSEA 集成包）MASLD vs MASH 髓系差异通路
# irGSEA 内部对 AUCell/ssGSEA/UGSEA/Seurat 打分做集成并给出差异统计
devtools::install_github("artem-sokolov/irGSEA")  # 实际仓库:  sishengjie/irGSEA 等，按原文引文确认
# 期望：angiogenesis、coagulation、hypoxia、炎症信号在 MASH 上调 (p<0.05~0.0001)
```

> irGSEA 为社区包，接口随版本变化；若版本不可得，用 GSVA(ssGSEA) + limma/wilcox 对 MASLD vs MASH 做通路差异，能得到等价结论（在文档中注明替代实现）。

---

## 6. 髓系拟时序与 State 图（Fig. 5A–D）

与 04 完全同构，只是对象换成髓系：

```r
# CytoTRACE：归一化 0–1；期望 VCAM1+/SAA1+ KC 高分（干性/祖态），DTNA+ Mac 低分
# Monocle2：根 = KC 富集 state；期望轨迹 VCAM1+ KC → DTNA+ Mac
# State 投影图（Fig.5B,C）：plot_cell_trajectory(color_by="State"/"subtype")
# 平滑热图（Fig.5D）：plot_pseudotime_heatmap 沿髓系分化程的基因动态
```

**核对（原文 Fig.5A–C）**：stemness 梯度 KC 高 → DTNA⁺ Mac 低；State 1→5 的演化方向与 KC→DTNA⁺ 一致。

---

## 7. 导师提示与陷阱

1. **DC 与巨噬的边界**：mDC/pDC 的 marker 表达稀疏，dropout 下打分易错；建议加 `pDC = LILRA4+IL3RA+CLEC4C` 三基因 AND 逻辑（至少 2/3 表达）复核。
2. **增殖群是状态不是身份**：MKI67⁺ Mac 与其他群可能重叠（双计数），若研究比例动态可保留，若做通讯分析建议保留但标注。
3. **M1/M2 是连续谱**：原文 Fig.4D 显示 DTNA⁺ Mac 同时有 M2 与促炎特征（"M1/M2 hybrid"），箱线图应如实呈现双向信号，不要强行二分类。
4. **DTNA 的表达稀疏性**：单细胞中 DTNA 检测率可能低，注释时以 cluster 平均表达 + ITGAX 共表达为准，不要用单细胞阈值筛"DTNA⁺ 细胞"做下游（会丢大量真阳）。
