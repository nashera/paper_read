# 09 bulk 转录组：ssGSEA 纤维化评分、PLG 高/低分层与乳酸代谢相关分析

> **前置依赖**：01（bulk 下载）。**对应原文**：Fig. 6E、Fig. 7A–C、Fig. S6、Fig. S7；Methods "Quantification of the fibrosis score"、"Statistical analysis"。
> **产出**：9 队列纤维化评分、PLG 75 分位分组的 DEG 与补体凝血 GSEA、乳酸代谢基因相关矩阵、小鼠队列 RUNX2 验证。

---

## 1. ssGSEA 原理

ssGSEA（Barbie 2009）是 GSEA 的单样本版本：对**单个样本**的基因表达做秩排序 $r_1 \ge r_2 \ge \cdots \ge r_N$（$N$ 为基因总数），计算某基因集 $S$（$|S| = N_H$）的富集分数。定义沿排序列表第 $j$ 位的命中/未命中累积量：

$$P_{\mathrm{hit}}(S, j) = \frac{\sum_{g \in S,\; \mathrm{rank}(g) \le j} |r_g|^{\alpha}}{\sum_{g \in S} |r_g|^{\alpha}}, \qquad P_{\mathrm{miss}}(S, j) = \frac{\sum_{g \notin S,\; \mathrm{rank}(g) \le j} 1}{N - N_H}$$

$$\mathrm{ES}(S) = \sum_{j=1}^{N} \left[ P_{\mathrm{hit}}(S, j) - P_{\mathrm{miss}}(S, j) \right]$$

其中 $\alpha \in \{0.75, 1\}$ 为加权指数（ssGSEA 常用 0.25 或 1，R 的 GSVA 默认 $\alpha$ 由参数 `tau` 控制）。ES 也可取累积曲线偏离 0 的最大绝对值（与 GSEA 一致的带符号版本）。

与 GSEA 的区别：ssGSEA 不需要表型分组，每个样本独立得分，因此可用于①连续评分（纤维化程度）②作为特征进入相关/回归模型。**原文用 MSigDB 的 HP_HEPATIC_FIBROSIS（Human Phenotype 分支）基因集对 bulk 样本打纤维化分。**

实现：R 包 GSVA 的 `ssgsea` 方法，或 Python `gseapy.ssgsea`。

---

## 2. 代码

### 2.1 队列表达矩阵标准化（含芯片探针映射）

```r
# code/09_bulk_ssgsea.R
library(GEOquery); library(GSVA); library(clusterProfiler)
nine <- c("GSE126848","GSE135251","GSE167523","GSE48452","GSE66676",
          "GSE89632","GSE164760","GSE162694","GSE24807")
expr <- list()
for (g in nine) {
  es <- getGEO(g, GSEMatrix = TRUE)[[1]]
  m  <- exprs(es)
  if (is_microarray(es)) {          # Affymetrix：探针→symbol，取每基因最大均值探针
    ann <- fvarLabels(es); sym <- fData(es)$`Gene Symbol` %||% annotate_probe(g)
    m <- aggregate_by_symbol(m, sym)
  }
  m <- log2(m + 1)                   # RNA-seq 若为 FPKM/TPM 同样 log2(x+1)
  expr[[g]] <- m
}
```

### 2.2 纤维化评分

```r
# MSigDB: HP_HEPATIC_FIBROSIS（https://www.gsea-msigdb.org/gsea/msigdb/cards/HP_HEPATIC_FIBROSIS）
hpc <- read.gmt("data/gmt/hp_hepatic_fibrosis.gmt")   # 或从 MSigDB 导出 h.all+hp 集合
fib_scores <- lapply(expr, function(m)
  ssGSEA(m, hpc$gene, method = "ssgsea", normalize = TRUE)$ssGSEA[1, ])
# 保存每队列样本级纤维化分
```

### 2.3 PLG 与纤维化/RUNX2 相关（与 07 文档共用）

```r
g <- "GSE126848"; m <- expr[[g]]
cor.test(fib_scores[[g]], m["PLG", ])      # 期望 R=0.639, p=1.1e-4
cor.test(m["RUNX2", ], m["PLG", ])         # 期望 R=0.655, p=6.33e-5
```

### 2.4 PLG 高/低分组（DTNA⁺ Mac 来源，75 分位）与 DEG（Fig. 7A）

原文逻辑：在**单细胞**的 DTNA⁺ Mac 中计算 PLG 表达，以 75 分位把**样本/细胞**分为 high-PLG 与 low-PLG，再比较 DTNA⁺ Mac 的 DEG。两种可复现口径（建议都做并报告）：

```r
mye <- readRDS("processed/myeloid_annotated.rds")
dtna <- subset(mye, subtype == "DTNA+ Mac")
thr  <- quantile(dtna[["PLG"]][,1], 0.75)          # 75th percentile（原文）
dtna$PLGgrp <- ifelse(dtna[["PLG"]][,1] > thr, "high", "low")
# 口径 A（细胞级 pseudo-bulk 更稳）：按 sample × PLGgrp 做 pseudo-bulk 后 DE
pb <- AggregateExpression(dtna, group.by = c("sample","PLGgrp"))$RNA
deg <- limma_voom_or_edgeR(pb_high_vs_low)         # high vs low
# 口径 B（细胞级 wilcox，原文 FindMarkers 风格）
deg <- FindMarkers(dtna, ident.1 = "high", ident.2 = "low", group.by = "PLGgrp")
fwrite(deg, "results/Fig7A_PLG_high_low_DEG.csv")
EnhancedVolcano(deg, lab = rownames(deg), x = "avg_log2FC", y = "p_val_adj")
```

### 2.5 补体与凝血级联 GSEA（Fig. 7B–C）

```r
rank <- sort(deg$avg_log2FC, decreasing = TRUE)
g <- gseKEGG(rank, organism = "hsa")               # 或 GSEA(HALLMARK/KEGG)
# 期望 top：Complement and coagulation cascades（hsa04610）在 high-PLG 上调
gseaplot2(g, geneSetID = which(g$Description == "Complement and coagulation cascades"))
```

### 2.6 乳酸代谢相关基因相关（Fig. S6）

```r
lac <- c("LDHA","LDHB","SLC16A1","SLC16A3","SLC2A1","HK2","PDK1","MCT1",
         "GPR81","HIF1A")   # 乳酸产生/转运/ sensing 基因集（按原文补充表核对）
key <- c("RUNX2","PLG","PARD3")
cm  <- cor(t(m[c(key, lac), ]), method = "pearson")[key, lac]
pheatmap(cm, display_numbers = TRUE)   # 期望三关键基因与多个乳酸基因显著相关
```

### 2.7 小鼠 HFHFD 队列 RUNX2 验证（Fig. S7A–B）

GSE166504（小鼠）中比较 30 周 HFHFD vs 对照的髓系 RUNX2 表达（wilcox），并在单细胞 UMAP/FeaturePlot 展示髓系 RUNX2 上调。

---

## 3. 陷阱

1. **芯片平台探针映射**：GSE48452/GSE66676 等 Affy 队列一个 symbol 对应多探针，取"该探针组中方差最大者"或"均值最大者"需固定一种并在代码注释；不同选择会小幅影响相关系数。
2. **ssGSEA 的 normalize 参数**：GSVA 的 `ssgsea.normalize=TRUE` 会把 ES 归一到 [0,1] 或排名尺度；跨队列比较时务必每队列内统一口径。
3. **75 分位的对象**：是"DTNA⁺ Mac 的 PLG 表达分布"的分位，作用于细胞或 pseudo-bulk 样本；不要误用全髓系或全肝的分位。
4. **相关≠因果**：本章全部为观察性相关；RUNX2→PLG 的因果证据来自 10 文档的 ATAC motif，报告中要分层表述。
