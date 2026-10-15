# 04 肝星状细胞（HSC）五亚型细分、差异分析与分化轨迹

> **前置依赖**：02。**对应原文**：Fig. 3 全图；Methods "Cell type annotation and subclustering"、"Developmental trajectory inference"；Suppl. Table 4–5。
> **产出**：5 亚型注释、Ro/e 动态、DEG 火山图 + KEGG、CytoTRACE 干性、Monocle2 拟时序轨迹、平滑热图。

---

## 1. 生物学背景

HSC 是肝纤维化的核心效应细胞：静息态（qHSC）储存维生素 A 脂滴；损伤后活化为 aHSC（肌成纤维细胞样），分泌 COL1A1/COL3A1 等 ECM。原文在单细胞层面进一步分辨出 **5 种状态**：

| 亚型 | 特征 marker（原文 Fig.3B / Suppl. Table 4） | 生物学含义 |
|---|---|---|
| qHSC | RELN 等静息相关基因 | 脂滴储存、静息 |
| inHSC | 中间态特征 | 活化过渡 |
| aHSC | KALRN、ITGBL1（炎症/活化） | 肌成纤维样、产 ECM |
| hep-like HSC | 肝细胞样表达程序 | 转分化/可塑性 |
| SMC | MYH11 等平滑肌分化基因 | 平滑肌样转化 |

---

## 2. 子聚类与注释

```r
# code/04_hsc_subtype.R
library(Seurat)
so <- readRDS("processed/integrated_harmony.rds")
hsc <- subset(so, celltype == "HSC")          # n ≈ 26,073

# 子集内重跑完整流程（关键：HVG 必须重算，全局 HVG 会漏掉亚型特异基因）
hsc <- NormalizeData(hsc) |> FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |> RunPCA(npcs = 30)
hsc <- RunHarmony(hsc, group.by.vars = "dataset", dims.use = 1:20)   # 子集仍含批次
hsc <- FindNeighbors(hsc, reduction = "harmony", dims = 1:20) |>
       FindClusters(resolution = 0.5) |>
       RunUMAP(reduction = "harmony", dims = 1:20, seed.use = 20260101)
DimPlot(hsc, label = TRUE)                     # 复现 Fig.3A

# --- 亚型 marker（原文 + 文献惯例）---
sub_mk <- list(
  qHSC        = c("RELN","LRAT","ADH1B","CILP"),          # RELN 为原文点名的静息基因
  inHSC       = c("PDGFRB","NOTCH3","ACTA2-low"),        # 过渡态：间质基因上调、ACTA2 中等
  aHSC        = c("KALRN","ITGBL1","ACTA2","TAGLN","COL1A1"),
  `hep-like`  = c("ALB","APOA2","TTR"),                  # 肝细胞样程序
  SMC         = c("MYH11","ACTG2","CNN1"))
# 用 02 文档的打分-最大法注释；cluster 数 ≠ 5 时合并/拆分至 marker 一致
DotPlot(hsc, features = unique(unlist(sub_mk))) + RotatedAxis()
ggsave("figures/Fig3B_hsc_dotplot.pdf")

# --- Ro/e（复用 03 的函数）：MASLD vs MASH；eMASH vs aMASH ---
# 期望结果：aHSC 的 Ro/e 随进展升高、qHSC 下降；aMASH 的 aHSC > eMASH（Fig.3C,D）
```

---

## 3. HSC 的 DEG 与通路（Fig. 3E–F）

```r
deg <- FindMarkers(hsc, ident.1 = "MASH", ident.2 = "MASLD",
                   group.by = "disease", test.use = "wilcox",
                   logfc.threshold = 0.25, min.pct = 0.1)
fwrite(deg, "results/Fig3E_HSC_DEG_MASH_vs_MASLD.csv")

# 火山图
library(EnhancedVolcano)
EnhancedVolcano(deg, lab = rownames(deg), x = "avg_log2FC",
                y = "p_val_adj", pCutoff = 0.05, FCcutoff = 0.5)

# KEGG 富集（上调基因）
up <- rownames(deg[deg$avg_log2FC > 0.25 & deg$p_val_adj < 0.05, ])
ego <- enrichKEGG(up, organism = "hsa", pAdjustMethod = "BH")
barplot(ego, showCategory = 15)
ggsave("figures/Fig3F_HSC_KEGG.pdf")
```

**核对**：原文报告 MASLD 组 HSC 富集 Axon Guidance、Relaxin Signaling、Dopaminergic Synapse；MASH 组富集 Cholesterol Metabolism、Bile Secretion。若方向/条目不一致，先检查 DEG 的 ident.1/ident.2 方向。

---

## 4. CytoTRACE：分化潜能的定量

### 4.1 原理

CytoTRACE（Gulati et al. 2020）的核心假设：**分化程度低的细胞表达更多种类的基因**（转录组复杂度高）。算法：

1. 对每个细胞 $i$ 计算表达基因数 $n_i = \#\{g : x_{gi} > 0\}$；
2. 对每个基因 $g$，计算其表达向量 $\mathbf{x}_g$ 与 $\mathbf{n}$ 的 Pearson 相关，得到基因-分化潜能关联：

$$\rho_g = \mathrm{cor}(\mathbf{x}_g, \mathbf{n})$$

3. 取 $\rho_g$ 最正相关的前 ~200 个基因组成集合 $G^*$，以其（秩标准化后）平均表达作为该细胞的 CytoTRACE 原始分：

$$C_i = \frac{1}{|G^*|}\sum_{g \in G^*} \tilde{x}_{gi}$$

4. 用马尔可夫平滑在 k-NN 图上迭代平滑（$C \leftarrow P\,C$，$P$ 为行归一化转移矩阵），减少 dropout 噪声；
5. 归一化到 $[0,1]$：

$$\hat{C}_i = \frac{C_i - \min_j C_j}{\max_j C_j - \min_j C_j}$$

$\hat{C}_i = 1$ 为最高分化潜能（最"干"），$\hat{C}_i = 0$ 为终末分化。

**为什么它能定轨迹起点**：Monocle2 的拟时序是无根树，方向需要外部信息；CytoTRACE 提供"哪端是祖细胞"的生物学锚点。原文正是用 CytoTRACE 确认 qHSC 位于高分端、aHSC/SMC 位于低分端。

### 4.2 代码

```r
library(cytoTRACE2)   # 或原版 CytoTRACE R 包
ct <- cytoTRACE2(as.matrix(GetAssayData(hsc, "counts")), platform = "10xv3")
hsc$cytotrace <- ct$ct_score          # 0–1
FeaturePlot(hsc, "cytotrace")        # 复现 Fig.3G：qHSC 高分、aHSC/SMC 低分
# 组间检验：wilcox 比较各亚型 cytotrace
```

---

## 5. Monocle2 拟时序

### 5.1 原理

Monocle2 流程：① 差异/高负载基因筛选 → ② **DDRTree**（Discriminative Dimensionality Reduction with Trees）降维：在低维空间同时学习一个**主图（principal graph，树结构）**与细胞到图的投影，目标函数为

$$\min_{\mathbf{B}, \mathbf{Y}, \mathbf{Z}} \; \left\| \mathbf{X} - \mathbf{B}\mathbf{Z} \right\|_F^2 \;+\; \gamma \left\| \mathbf{Y} - \mathbf{Z} \right\|_F^2 \;+\; \lambda \sum_{(j,k) \in \mathcal{E}} \left\| \mathbf{y}_j - \mathbf{y}_k \right\|^2$$

其中 $\mathbf{B}$ 为投影基、$\mathbf{Z}$ 为细胞在主图节点上的软分配、$\mathbf{Y}$ 为主图节点坐标、$\mathcal{E}$ 为树的边集（树脊惩罚保证主图是树而非任意图形）；③ 沿主图测地距离排序得 pseudotime $t_i = d_{\mathrm{geodesic}}(\mathbf{z}_i, \mathbf{z}_{\mathrm{root}})$；④ 以指定根（root）定向。

与 PCA 的区别：PCA 只能得到线性子空间；DDRTree 得到的是**弯曲的树状流形**，能表达分叉（branching）——这正是"qHSC → inHSC → 三个终末分支"的结构。

### 5.2 代码

```r
library(monocle)
cds <- newCellDataSet(as(as.matrix(GetAssayData(hsc, "counts")), "dgCMatrix"),
                      phenoData = new("AnnotatedDataFrame", data = hsc@meta.data),
                      featureData = new("AnnotatedDataFrame",
                                        data = data.frame(gene_short_name = rownames(hsc),
                                                          row.names = rownames(hsc))),
                      lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds) |> estimateDispersions()

# 用高变/差异基因做排序基因（Monocle2 惯例）
dg <- differentialGeneTest(cds, fullModelFormulaStr = "~subtype_label")
ord  <- row.names(subset(dg, qval < 0.01))
cds  <- setOrderingFilter(cds, ord)
cds  <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds  <- orderCells(cds)                       # 无根伪时序
# 定向：以 CytoTRACE 最高的细胞（qHSC 富集区）为根
root <- rownames(pData(cds))[which.max(hsc$cytotrace[match(rownames(pData(cds)), colnames(hsc))])]
cds  <- orderCells(cds, root_state = cds@auxOrderData$HS$root_state_from_cell(root))

plot_cell_trajectory(cds, color_by = "subtype_label")   # Fig.3H
plot_cell_trajectory(cds, color_by = "State")           # Fig.3I：State 4/5=qHSC 祖端, State 1=终末
plot_pseudotime_heatmap(cds, num_clusters = 4, show_rownames = FALSE)  # Fig.3J 平滑热图
```

**核对点**：
- 轨迹拓扑 = qHSC → inHSC → {aHSC, hep-like HSC, SMC} 三分叉；
- State 编号语义：原文 State 4/5 对应 qHSC 祖态、State 1 对应 aHSC/SMC 终末（Monocle 的 State 编号是自动的，需按亚型占比人工对应，不要硬编码编号）；
- 平滑热图（plot_pseudotime_heatmap 内部按伪时序做平滑聚类）呈现沿程动态基因模块。

### 5.3 把 CytoTRACE 叠到 Monocle 对象（原文做法）

```r
pData(cds)$cytotrace <- hsc$cytotrace[match(rownames(pData(cds)), colnames(hsc))]
plot_cell_trajectory(cds, color_by = "cytotrace")
# 检查：pseudotime 与 cytotrace 应负相关（Spearman），即越"干"越早
```

---

## 6. 导师提示

1. **Monocle2 vs Monocle3**：原文指定 Monocle2（DDRTree）。Monocle3 用 UMAP+Simple Principal Graphs，拓扑更灵活但不可直接互换；复现请忠于原文。
2. **排序基因选择决定轨迹**：用全部基因会被管家基因主导；用亚型 DEG（上文 differentialGeneTest）是标准做法。
3. **根的生物学论证**：不能只靠 CytoTRACE 一个证据；结合 qHSC 的静息 marker（RELN/LRAT）在根端富集做交叉验证。
4. **分支点基因**：可用 `BEAM(cds, branch_point=..., branch_states=c(a,b))` 提取决定 aHSC vs SMC 命运的分支依赖基因，作为机制候选（原文未展开，是好的扩展点）。
