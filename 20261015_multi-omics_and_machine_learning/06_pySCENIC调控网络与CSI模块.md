# 06 pySCENIC 转录调控网络推断与 CSI regulon 模块

> **前置依赖**：05（髓系注释对象）。**对应原文**：Fig. 5E–G、Fig. S3；Methods "Transcription factor regulatory network analysis"、"Regulon module identification"；Suppl. Table 7。
> **产出**：360 个显著 regulon、M1–M17 模块、DTNA⁺ Mac 富集模块（M2/M3/M10/M12，74 TF）、RSS Top10、Venn 交集 {BHLHE41, RUNX2, BACH2, MEF2C}。

---

## 1. SCENIC 三步法的原理

SCENIC（Aibar 2017）从表达数据重建 **regulon = 转录因子 TF + 其直接靶基因集合**，并在单细胞水平量化每个 regulon 的活性。三步：

### 步骤 1：共表达推断（GENIE3 / GRNBoost2）

把"预测每个基因 $g$ 的表达"建模为回归/集成树问题：以所有 TF 的表达 $\{x_{t,i}\}_{t \in \mathrm{TF}}$ 为特征，预测靶基因 $g$ 的表达 $x_{g,i}$。GENIE3 用 **随机森林/ExtraTrees** 对每个靶基因拟合：

$$\hat{x}_g = F_g\left(\mathbf{x}_{\mathrm{TF}}\right) = \sum_{t=1}^{T} \alpha_t\, f_t\left(\mathbf{x}_{\mathrm{TF}}\right)$$

权重 $w_{tg}$（TF $t$ → 靶基因 $g$ 的调控强度）取特征重要性，即 $t$ 在全部树中造成的**不纯度下降总和**：

$$w_{tg} = \sum_{\text{树}} \sum_{\text{以 } t \text{ 分裂的节点}} \left( n_{\text{node}} \cdot \Delta \mathrm{impurity} \right)$$

GRNBoost2 是其梯度提升（GBM）加速版，数学目标等价：

$$\min_{F_g} \; \sum_i \left( x_{g,i} - F_g(\mathbf{x}_{\mathrm{TF},i}) \right)^2 \;+\; \Omega(F_g)$$

输出：TF×靶基因 的权重矩阵（co-expression module）。

> 关键局限：共表达 ≠ 直接调控（可能经中间因子）。因此需要步骤 2 的 motif 证据过滤。

### 步骤 2：motif 富集剪枝（RcisTarget）

对每个 TF 的候选靶集 $T$，检查其 **顺式调控序列（启动子 ±、或增强子区域）中该 TF 的结合 motif 是否富集**：
- 预计算全基因组（人 hg38）每个 motif $m$ 在每条序列上的排名（ROC 曲线下面积 AUC 排序数据库，如 `hg38-500bp-upstream-10x10species` / `hg38-tss-centered-10kb`）。记 $N$ 条参考序列按 motif $m$ 的打分降序排名 $1 \dots N$，$N_r^m$ = 前 $r$ 名中属于候选靶集 $T$ 的序列数，富集分数为恢复曲线 AUC 的归一化偏差：

$$\mathrm{NES}(m, T) = \frac{\frac{1}{N}\sum_{r=1}^{N} N_r^m \;-\; \mu}{\sigma}, \qquad \mu = \frac{|T|(N+1)}{2}, \;\; \sigma = \sqrt{\frac{|T|(N-|T|)(N+1)}{12}}$$

- NES 显著（FDR 校正）者保留该 motif；
- 只保留"motif 富集显著"的 TF-靶边 → 得到 **regulon**（直接靶版本 direct，+ 扩展的间接靶版本标记为 extended）。

### 步骤 3：单细胞 regulon 活性（AUCell）

对每个细胞 $i$，把基因按表达降序排名，regulon $r$ 的活性为靶基因在排序前列 $[1, T_{\max}]$ 的恢复曲线面积：

$$\mathrm{AUCell}(r, i) = \frac{\sum_{j=1}^{T_{\max}} N_r^m(j) \cdot \mathbb{1}[g_j \in r]}{\text{归一化项}} \in [0, 1]$$

（实际实现为 $\sum_{j} \frac{1}{|r \cap L_{\le j}|}$ 型归一，等价于"regulon 基因在排名前列的富集面积"。）AUC 高 = 该 regulon 在此细胞"开启"。对每个 regulon 设 AUC 阈值二值化 → regulon 活性矩阵（细胞 × regulon，0/1 或连续 AUC）。

---

## 2. pySCENIC 执行（Python，CLI 三步）

```bash
# code/06_scenic.sh
# 输入：raw count 矩阵（原文：raw count matrix as input）+ 细胞注释
# 准备 expr.mat.gz（基因×细胞，raw counts）与 mm 数据库
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38-500bp-upstream-10x10species.mc_v10_clust.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# 步骤1+2：GRNBoost2 共表达 + RcisTarget motif 剪枝 → regulon(.adjacencies)
pyscenic grnboost expr.mat.gz \
    -o adjacencies.csv --num_workers 32 \
    --expression_mtx_header gene cell   # 按实际矩阵方向调整
# （grnboost 等价命令: pyscenic grn）

pyscenic ctx adjacencies.csv \
    hg38-500bp-upstream-10x10species.mc_v10_clust.genes_vs_motifs.rankings.feather \
    hg38-tss-centered-10kb-10x10species.mc_v10_clust.genes_vs_motifs.rankings.feather \
    --motif_annotations motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --num_workers 32 -o regulons.csv \
    --mask_dropouts

# 步骤3：AUCell 单细胞活性
pyscenic aucell expr.mat.gz regulons.csv -o aucell_mtx.csv --num_workers 32
```

**规模提示**：3.8 万髓系细胞 × 全 TF 的 GRNBoost 在 32 核约 1–2 天；可先对髓系子集运行（原文 regulon 分析针对髓系/各细胞群）。`--mask_dropouts` 按原文流程保留 dropout 掩码处理。

---

## 3. 显著 regulon 筛选与 CSI 模块聚类

### 3.1 CSI（Connection Specificity Index）原理

CSI（Suo et al. 2018, Cell Rep）用于度量"两个 regulon 的活性模式是否在上下文（所有其他 regulon）中特异性地相似"：

1. 计算 regulon 活性矩阵（细胞 × regulon）两两间的 **Pearson 相关** $\mathrm{PCC}(A,B)$；
2. 对固定对 $(A,B)$：

$$\mathrm{CSI}(A,B) = \frac{1}{n-2}\left|\left\{ C \neq A,B \;:\; \mathrm{PCC}(A,C) < \mathrm{PCC}(A,B) \;\wedge\; \mathrm{PCC}(B,C) < \mathrm{PCC}(A,B) \right\}\right|$$

即"A 与 B 的相关，比它们各自与几乎所有其他 regulon 的相关都高"的 regulon 占比。$\mathrm{CSI} \in [0,1]$，高 CSI = 特异性伙伴。**优点**：普通相关会被"全局共变模块"污染（所有活跃 regulon 互相关高），CSI 用相对排名消除这种 hub 效应。

3. 对 CSI 矩阵以 **欧氏距离做层次聚类**（原文）→ regulon 模块 M1–M17：

$$d(\mathbf{CSI}_A, \mathbf{CSI}_B) = \sqrt{\sum_{j}\left(\mathrm{CSI}(A,j) - \mathrm{CSI}(B,j)\right)^2}$$

4. 以 $\mathrm{CSI} > 0.7$ 为边阈值构建 regulon 关联网络可视化；
5. 模块在细胞类型 $k$ 的活性 = 模块内成员 regulon 在该类型所有细胞上的平均 AUCell 活性：

$$\mathrm{Act}(M, k) = \frac{1}{|M|}\sum_{r \in M} \frac{1}{|C_k|}\sum_{i \in C_k} \mathrm{AUCell}(r, i)$$

取每模块活性最高的细胞类型（top-ranking）。

### 3.2 代码（R 或 Python 均可，此处 R）

```r
# code/06_csi_modules.R
auc <- as.matrix(read.csv("aucell_mtx.csv", row.names = 1))   # 细胞 × regulon
auc <- auc[, colSums(auc > 0) > 0.01 * nrow(auc)]             # 显著/活跃 regulon 过滤
# 原文最终 360 个显著 regulon：按"在≥1% 细胞活跃 + AUCell 阈值二值化后非空"筛选，
# 并记录筛选前后的 regulon 数以复现 360 这个数字

pcc <- cor(auc, method = "pearson")                           # regulon × regulon
n   <- ncol(auc)
csi <- matrix(0, n, n, dimnames = dimnames(pcc))
for (i in 1:(n-1)) for (j in (i+1):n) {
  a <- pcc[i, j]
  ok <- (pcc[i, -c(i,j)] < a) & (pcc[j, -c(i,j)] < a)
  csi[i, j] <- csi[j, i] <- mean(ok)
}
diag(csi) <- 1

# 层次聚类（欧氏距离）→ 切 17 模块
hc  <- hclust(dist(csi, method = "euclidean"), method = "average")
mod <- cutree(hc, k = 17)                                      # 原文 M1–M17
table(mod)

# CSI>0.7 关联网络
library(igraph)
el <- which(csi > 0.7 & upper.tri(csi), arr.ind = TRUE)
g  <- graph_from_data_frame(data.frame(from = rownames(csi)[el[,1]],
                                       to   = colnames(csi)[el[,2]]), directed = FALSE)
plot(g, vertex.size = 4)

# 模块 × 细胞类型活性热图（Fig.5E）
meta <- readRDS("processed/myeloid_annotated.rds")@meta.data
act  <- sapply(split(colnames(auc), mod), function(r)
          tapply(colMeans(auc[meta$cell_barcode, r, drop=FALSE]), meta$subtype, mean))
pmat <- t(scale(t(act)))
pheatmap(pmat, color = colorRampPalette(c("blue","white","red"))(100))
ggsave("figures/Fig5E_regulon_module_heatmap.pdf")
# 核对：DTNA+ Mac 富集模块 = M2, M3, M10, M12；四模块 TF 并集 = 74
tf_of <- function(m) unique(sub("\\(\\+\\)$|\\(-\\)$", "", names(mod)[mod == m]))
length(unique(unlist(lapply(c(2,3,10,12), tf_of))))   # 期望 74
```

---

## 4. RSS 特异性评分与 Top10 TF（Fig. 5F）

**RSS（Regulon Specificity Score）**：把某细胞类型 $k$ 的指示向量（$k$ 内细胞=1，其余=0）与每个 regulon 的活性向量比较，用**信息论的 Jensen–Shannon 散度**度量特异性（Su et al. 2016 的原始定义针对基因；SCENIC 社区实现将其用于 regulon 活性）：

$$\mathrm{RSS}(r, k) = 1 - \mathrm{JSD}(\mathbf{p}_r \,\|\, \mathbf{e}_k) = 1 - \frac{1}{2}D_{\mathrm{KL}}(\mathbf{p}_r \|\mathbf{m}) - \frac{1}{2}D_{\mathrm{KL}}(\mathbf{e}_k \|\mathbf{m}), \quad \mathbf{m} = \frac{\mathbf{p}_r + \mathbf{e}_k}{2}$$

其中 $\mathbf{p}_r$ 为 regulon $r$ 的活性在各细胞上的归一化分布（$\sum_i p_{ri}=1$），$\mathbf{e}_k$ 为类型 $k$ 的 0/1 指示归一化分布。$\mathrm{RSS} \to 1$ 表示 $r$ 的活性完全集中于类型 $k$。

```r
# 对每个髓系亚群计算全部 regulon 的 RSS，取 Top10 画 bar/dot（Fig.5F, S3B）
rss <- function(act_vec, ind) {
  p <- act_vec / sum(act_vec); q <- ind / sum(ind)
  m <- (p + q) / 2
  jsd <- 0.5*sum(p*log(p/m)) + 0.5*sum(q*log(q/m))
  1 - jsd
}
```

---

## 5. Venn 交集定位关键 TF（Fig. 5G）

两个 regulon 集合：
- 集合 1 = DTNA⁺ Mac 富集模块（M2/M3/M10/M12）的成员 TF；
- 集合 2 = DTNA⁺ Mac 的 RSS Top（高特异性）TF；

交集 = **BHLHE41(+), RUNX2(+), BACH2(+), MEF2C(+)**（"+"表示 regulon 为激活型/扩展型标记）。

```r
library(ggVennDiagram)
ggVennDiagram(list(Module_TFs = set1, RSS_TFs = set2), label_alpha = 0)
# 核对交集恰为上述 4 个 TF
```

**机制衔接（原文）**：RUNX2 已知通过 ITGAV 促 HSC 活化与纤维化、上调 MCP-1 促巨噬浸润、驱动肝细胞死亡——因此选 RUNX2 作为轴的上游 TF 进入 07/10 文档验证。

---

## 6. 导师提示

1. **360 这个数字的复现敏感性**：regulon 数取决于①活跃阈值②motif 数据库版本③输入细胞集。复现时记录三个版本参数；若得 300–400 属正常波动，但模块结构（DTNA⁺ Mac 富集模块）应稳定。
2. **GRNBoost vs GENIE3**：原文 Methods 写 GENIE3（概念）+ GRNboost（实现），二者目标函数同源；不要混用 pyscenic 的 `grn` (GENIE3) 与 `grnboost` 而不记录。
3. **CSI 的计算复杂度** O(n³)（n=regulon 数），360² 对级别可接受；向量化实现避免三重循环。
4. **regulon 命名**：pySCENIC 输出形如 `RUNX2(+)"或"RUNX2(extended)"`；统计 TF 数前必须剥离后缀，否则 74 会算错。
