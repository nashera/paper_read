# 07 CellChat 细胞通讯分析与 RUNX2–PLG–PARD3 轴验证

> **前置依赖**：04（HSC 亚型）、05（髓系亚群）、06（RUNX2）。**对应原文**：Fig. 6A–F、Fig. S4；Methods "Cell-cell communication analysis"。
> **产出**：MASLD/MASH 通讯网络对比、差异互作热图、PLG–PARD3 配受体证据、bulk 相关（R=0.639 / 0.655）、9 队列 DTNA⁺Mac–aHSC 相关矩阵。

---

## 1. CellChat 的数学原理

CellChat（Jin 2021/2025 Nat Protoc）从单细胞表达推断**配体–受体（L–R）介导的通讯概率**：

1. **输入**：归一化表达矩阵 + 细胞群标签；内置数据库 CellChatDB（人：~1500+ L–R 对，含分泌信号、ECM–受体、细胞接触三类）。
2. **通讯概率模型**：对配体 $L$（发送群 $i$）与受体 $R$（接收群 $j$），CellChat 定义通讯概率

$$P_{ij}^{(L,R)} \propto \left( N_L^{(i)} \cdot \bar{E}_L^{(i)} \right) \times \left( N_R^{(j)} \cdot \bar{E}_R^{(j)} \right)$$

其中 $N_L^{(i)}$ 为群 $i$ 中表达 $L$ 的细胞占比，$\bar{E}_L^{(i)}$ 为其平均表达量（受体侧同理）。实际实现中 $\bar{E}$ 用 **截尾均值**（truncatedMean，去掉最高/最低各 10%）以降低离群细胞影响。多亚基受体 $R = R_1 \oplus \cdots \oplus R_n$ 按数据库注释取乘积或均值：

$$\bar{E}_R^{(j)} = \prod_{m=1}^{n} \bar{E}_{R_m}^{(j)} \quad \text{或} \quad \frac{1}{n}\sum_{m=1}^{n} \bar{E}_{R_m}^{(j)}$$

并引入**辅因子修正**：共刺激受体/辅受体项相乘、抑制因子项相除，

$$P_{ij}^{\mathrm{final}} = P_{ij}^{(L,R)} \times \prod_{a \in \mathrm{co\text{-}activator}} \bar{E}_a^{(j)} \Big/ \prod_{b \in \mathrm{inhibitor}} \left(1 + \bar{E}_b^{(j)}\right)$$

3. **统计显著性**：置换检验——随机打乱细胞群标签重算概率，重复 $B$ 次得到经验 $p$ 值 $p = (1 + \#\{b: P^{*(b)} \ge P\})/(1+B)$，过滤随机接触。
4. **通路级汇总**：把 L–R 对映射到信号通路（如 TGFb、COLLAGEN、COMPLEMENT），通路级通讯强度 = 该通路下所有 L–R 对概率之和 $I_{ij}^{(w)} = \sum_{(L,R) \in w} P_{ij}^{(L,R)}$。
5. **组间比较**：对两个 CellChat 对象（MASLD vs MASH）做 `computeNetSimilarity`/差异分析，输出互作数与强度的差值热图 $\Delta_{ij} = I_{ij}^{\mathrm{MASH}} - I_{ij}^{\mathrm{MASLD}}$（原文 Fig.6B，范围 $-20 \sim 20$）。

---

## 2. 代码（R）

```r
# code/07_cellchat.R
library(CellChat)
so  <- readRDS("processed/integrated_harmony.rds")
# 通讯分析需要"发送/接收群"同时包含髓系亚群与 HSC 亚群：
# 用 04/05 的亚型标签替换七大类中的 MonoMac 与 HSC
meta <- so@meta.data
meta$group <- ifelse(meta$celltype == "MonoMac", meta$myeloid_sub,
              ifelse(meta$celltype == "HSC", meta$hsc_sub, as.character(meta$celltype)))
so$group <- meta$group

build <- function(dz) {
  idx <- so$disease == dz
  data <- GetAssayData(so, "data")[, idx]
  cc <- createCellChat(object = data, meta = so@meta.data[idx, ], group.by = "group")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)                       # 只保留数据库基因，提速
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}
cc_masld <- build("MASLD"); cc_mash <- build("MASH")

# Fig.6A：互作数网络图（圈大小=细胞数/连线=互作数）
par(mfrow = c(1,2))
netVisual_circle(cc_masld@net$count, vertex.weight = table(cc_masld@idents))
netVisual_circle(cc_mash@net$count,  vertex.weight = table(cc_mash@idents))

# Fig.6B：差异互作热图（MASH - MASLD），计数与强度两版
net_diff <- cc_mash@net$count - cc_masld@net$count
pheatmap(net_diff, color = colorRampPalette(c("blue","white","red"))(100),
         breaks = seq(-20, 20, length.out = 101))
ggsave("figures/Fig6B_diff_interaction.pdf")
# 核对：DTNA+ Mac 作为 sender → aHSC/hep-like/inHSC 的边在 MASH 为正差

# Fig.6C：DTNA+ Mac 发出的 L-R circle plot
netVisual_aggregate(cc_mash, signaling = "PLG", vertex.receiver = NULL)
# 或直接提取 L-R 级：
lr <- subset(cc_mash@net$prob, ...)  # 用 netVisual_circle(cc_mash@net$prob, ...) 全局
# 关注 PLG–PARD3 对在 DTNA+ Mac → HSC 方向的概率

# Fig.6D：PLG / PARD3 表达小提琴（按群 × 疾病）
VlnPlot(so, features = c("PLG","PARD3"), group.by = "group", split.by = "disease")
# 期望：PLG 在 MASH 的 DTNA+ Mac 激活；PARD3 在两种条件下均于 HSC 表达
```

---

## 3. bulk 验证：纤维化评分–PLG、RUNX2–PLG 相关（Fig. 6E）

纤维化评分用 ssGSEA（详见 09 文档），在 GSE126848 上：

```r
# code/07_bulk_cor.R（依赖 09 的 ssGSEA 函数）
fib <- ssgsea_score(expr_GSE126848, geneset = "HP_HEPATIC_FIBROSIS")
cor.test(fib, expr["PLG", ],  method = "pearson")   # 期望 R=0.639, p=1.1e-4
cor.test(expr["RUNX2", ], expr["PLG", ], method = "pearson")  # 期望 R=0.655, p=6.33e-5
# 散点 + 线性拟合 + 置信带，复现 Fig.6E
```

> 正态性检查：原文统计约定为"正态→Pearson，否则 Spearman"；对每个相关先做 Shapiro–Wilk，并在结果表中记录所用方法。若复现值与原文偏差 >0.1，检查①基因探针映射②样本过滤③log 转换口径。

**Fig. S4**：在其余 bulk 队列（GSE135251 等）重复上述两个相关，作为跨队列稳健性。

---

## 4. 九队列 DTNA⁺Mac–aHSC 相关矩阵（Fig. 6F）

原理：用单细胞参考的亚群 signature（DTNA⁺ Mac 与 aHSC 的 top marker 均值或 ssGSEA 分）对 bulk 样本打分，再在每队列内计算两 score 的 Pearson/Spearman 相关。

```r
sig_mac <- top_markers_myeloid[["DTNA+ Mac"]]   # 05 的 FindAllMarkers 输出 top50
sig_ahsc <- top_markers_hsc[["aHSC"]]
for (gse in nine_bulk) {
  s1 <- ssgsea_score(expr[[gse]], sig_mac)
  s2 <- ssgsea_score(expr[[gse]], sig_ahsc)
  r  <- cor.test(s1, s2)   # 期望 9 队列中 6 个 p<0.05 且 R>0.8（原文口径）
}
# 汇总为相关矩阵/热图（Fig.6F）；讨论中总相关 R=0.36, p<2.2e-16 为合并口径
```

**口径提醒**：原文摘要级数字有两处——Fig.6F 为"6/9 数据集 p<0.05, R>0.8"，Discussion 为合并相关 R=0.36。复现时分别记录两口径，不要混用。

---

## 5. 导师提示

1. **truncatedMean vs mean**：CellChat 默认三均值（去掉最高最低各 10%）稳健估计群内表达；换 mean 会放大高表达离群细胞的影响。复现用默认。
2. **min.cells=10** 过滤稀疏群伪互作；髓系亚群细胞数充足，不需放宽。
3. **PLG–PARD3 的生物学谨慎性**：PLG（纤溶酶原）经典为肝细胞分泌；原文主张 MASH 中 DTNA⁺ Mac 来源的 PLG 异常激活。复现者应同时报告肝细胞的 PLG 表达作为对照（ Violin 中肝细胞应高表达 PLG），并在讨论中承认该轴为"计算推断 + ATAC 佐证"，缺蛋白水平验证（原文自承局限）。
4. **方向性**：CellChat 的 sender/receiver 由 L/R 表达决定；PLG 为配体在 Mac、PARD3 为受体在 HSC，与原文 violin 一致才可宣称轴方向。
