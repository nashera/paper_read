# 08 空间转录组 RCTD 反卷积与巨噬细胞–aHSC 空间共定位

> **前置依赖**：02/05（参考图谱）、01（Visium 数据）。**对应原文**：Fig. 2I、Fig. 7D–E、Fig. S2、S9、S10；Methods "Robust cell type decomposition"。
> **产出**：7 张 MASH Visium 的 spot 级细胞类型比例、mono/mac–HSC 共定位图、DTNA⁺Mac–aHSC signature 相关散点。

---

## 1. RCTD 的统计原理

Visium 每个 spot（55 μm）含 1–10 个细胞，表达是**多种细胞类型的混合**。RCTD（Cable 2022, Nat Biotechnol）用单细胞参考做**带技术噪声建模的反卷积**：

1. **参考构建**：以 scRNA-seq（原文用验证队列 GSE166504 的 38 样本）计算每类细胞 $c$ 的平均表达谱 $\mu_{cg}$ 与**基因特异性 overdispersion**（技术方差模型）；
2. **spot 模型**：spot $s$ 对基因 $g$ 的计数服从混合泊松：

$$y_{sg} \sim \mathrm{Poisson}\!\left(\mathrm{lib}_s \sum_{c} w_{sc}\,\mu_{cg}\right), \qquad \sum_c w_{sc} = 1,\; w_{sc} \ge 0$$

其中 $\mathrm{lib}_s$ 为 spot 文库大小，$w_{sc}$ 为细胞类型权重。RCTD 的关键是在似然中同时估计**平台特异性噪声参数**（基于参考单细胞的基因级 overdispersion 与 UMI 深度校正项），从而在 Visium（低深度、高 dropout）与 10x 单细胞（参考）之间做平台桥接，用 EM 类迭代求解 $w$；
3. **模式**：`doublet`（每个 spot 最多 2 类，即 $w_{sc}$ 至多 2 个非零元）与 `full`（多类 + 稀疏惩罚）。原文对全部 Visium 做 spot 注释 + 单细胞分辨率 deconvolution，建议 full mode；
4. **输出**：spot × celltype 的权重矩阵（细胞数估计 $\hat{n}_{sc} = w_{sc} \times \hat{N}_s$，$\hat{N}_s$ 为 spot $s$ 的细胞总数估计）。

**共定位定量**：对每个 spot 得到 w(DTNA⁺Mac/mono-mac) 与 w(aHSC/HSC)，计算两 signature score 的 Pearson 相关 + 空间自相关对照（打乱 spot 坐标的零分布），即原文 Fig.7D–E 的散点与 p<0.001。

---

## 2. 代码（R，spacexr）

```r
# code/08_rctd.R
library(spacexr)
# --- 参考图谱：GSE166504（原文指定），按 02 流程注释到"主要类型(+亚群)"级 ---
ref <- readRDS("processed/ref_GSE166504_annotated.rds")
counts_ref  <- as.matrix(GetAssayData(ref, "counts"))
cell_types  <- ref$celltype           # 或更细的亚群标签
nUMI_ref    <- colSums(counts_ref)
reference   <- Reference(counts_ref, cell_types, nUMI_ref)

# --- Visium spot 数据（GSE248077，7 个 MASH 样本）---
samples <- list.dirs("raw/spatial/GSE248077", recursive = FALSE)
for (s in samples) {
  mtx  <- Read10X_spatial(file.path(s, "filtered_feature_bc_matrix"))
  coord <- read.csv(file.path(s, "tissue_positions_list.csv"), header = FALSE)
  coords <- as.matrix(coord[, c(3,4)])       # row/col 空间坐标（按实际列序）
  nUMI_s <- colSums(mtx)
  spot   <- SpatialRNA(coords, mtx, nUMI_s)

  rctd <- create.RCTD(spot, reference, max_cores = 8, CELL_MIN = 200)
  rctd <- run.RCTD(rctd, doublet_mode = "full")

  w <- rctd@results$weights                # spot × celltype 稀疏权重
  fwrite(as.data.frame(as.matrix(w)), paste0("results/rctd_weights_", basename(s), ".csv"))
}
```

### 2.1 共定位与 signature score（Fig. 2I / 7D–E）

```r
# 方法 A（权重级）：spot 的 mono/mac 权重 vs HSC 权重
cor.test(w[,"MonoMac"], w[,"HSC"], method = "pearson")

# 方法 B（signature score 级，原文 Fig.2I/7D 的 "signature score" 图）：
# 用 Seurat AddModuleScore 在 Visium 表达上直接打分
vis <- Read10X_spatial(...) |> CreateSeuratObject()
vis <- AddModuleScore(vis, features = list(mac = sig_mac, hsc = sig_ahsc))
SpatialFeaturePlot / FeaturePlot(vis, c("mac1","hsc1"))     # Fig.2I 风格并列图
ggplot(scores, aes(mac1, hsc1)) + geom_point() + geom_smooth(method = "lm") +
  stat_cor()   # 期望显著正相关 p<0.001（Fig.7E）

# 空间零模型：置换 spot 坐标 1000 次重算相关，得到经验 p（比解析 p 更稳）
```

### 2.2 PLG–PARD3 空间共表达（Fig. S10）

```r
# 在 Visium 上计算 PLG 与 PARD3 的 spot 级相关 + 共表达热点图
cor.test(vis[["PLG"]], vis[["PARD3"]])
SpatialDimPlot 叠加热图展示共表达区域
```

---

## 3. 核对与陷阱

1. **参考物种一致性**：GSE166504 为小鼠（HFHFD），Visium GSE248077 为人 MASH——原文确以小鼠队列作参考（Methods 明写）。跨物种反卷积依赖直系同源基因保守性，**必须做 ortholog 映射后仅保留 1:1 直系基因**，并在复现报告中声明该局限；若可得人肝 Visium 参考，作为敏感性分析补一组。
2. **CELL_MIN / max_cores**：参考每类细胞 ≥200 个才稳定；髓系亚群在参考中若过稀，合并到 mono/mac 级再做亚群共定位（用 signature score 法 B 补足亚群分辨率）。
3. **full vs doublet**：full mode 允许 >2 类，肝组织 spot 常含肝细胞+内皮+免疫，必须 full。
4. **坐标方向**：Visium 的 row/col 与图像 y 轴可能翻转；叠加组织图时核对形态（肝小叶/门管区位置）。
