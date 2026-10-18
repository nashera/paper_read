# 04 精度加权 LMM 与经验贝叶斯 moderated t

> **目标文献**：Hoffman GE, et al. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* Nature Communications (2026) 17:9597. DOI: 10.1038/s41467-026-75680-8
> **模块对应原文**：Methods「Precision-weighted linear mixed models」「Empirical Bayes shrinkage for precision-weighted linear mixed models」；Fig. 1；引文 23（dream, Hoffman & Roussos 2021）、24（voom）、25（Smyth 2004）、30（lme4）
> **前置依赖**：`03_两阶段精度权重与均值方差趋势.md`（`processed/proc_<dataset>.rds`）
> **复现难度 / 预计工时 / 计算资源**：★★★★★ / 3–5 人日（含推导与单元测试）/ 内存 32–64 GB；全队列 DE 数小时

---

## ① 科学问题与逻辑位置

- **上游**：`dreamletProcessedData`（log2 CPM `E` + 权重 `weights` + 样本元数据）。
- **下游**：`dreamletResult` / `MArrayLM2`——包含每个基因的 $\hat\beta$、标准误、t、p、$d_0$、$s_0^2$、残差自由度。05 号抽取结果，06 号复用同一模型做方差分解，07 号用 t 统计量做基因集检验。
- **本模块回答三个问题**：
  1. **模型**：为什么用线性混合模型（LMM）而不是固定效应、不是负二项 GLMM？
  2. **权重**：精度权重怎样进入 LMM（加权最小二乘）？
  3. **小样本**：只有几十个受试者时，逐基因方差估计很不稳；经验贝叶斯怎样"借基因间的信息"把它稳住，并把零分布修正成 moderated t？

---

## ② 背景原理

### 2.1 线性混合模型（LMM）的设定

对基因 $g$、样本 $s$：

$$
y_{g,s} = \mathbf{x}_s^\top \boldsymbol\beta_g + \mathbf{z}_s^\top \mathbf{b}_g + \varepsilon_{g,s},
\qquad
\mathbf{b}_g \sim \mathcal N(\mathbf 0, \mathbf D_g),
\qquad
\varepsilon_{g,s}\sim\mathcal N\!\left(0,\frac{\sigma_g^2}{w_{g,s}}\right).
$$

- $\mathbf x_s$：固定效应设计（疾病状态、年龄、性别……），$\boldsymbol\beta_g$ 是我们**要检验**的。
- $\mathbf z_s$：随机效应设计（受试者、pool、批次……），$\mathbf b_g$ 是**要吸收**的。
- $w_{g,s}$：03 号的精度权重；方差与权重成反比，故高精度观测在拟合中占更大权重。

写成矩阵：

$$
\mathbf y_g = \mathbf X\boldsymbol\beta_g + \mathbf Z\mathbf b_g + \boldsymbol\varepsilon_g,\qquad
\operatorname{Var}(\mathbf y_g)=\mathbf V_g=\sigma_g^2\left(\mathbf W_g^{-1}+\mathbf Z\mathbf D_g^{*}\mathbf Z^\top\right),
$$

其中 $\mathbf W_g=\operatorname{diag}(w_{g,s})$，$\mathbf D_g^{*}=\mathbf D_g/\sigma_g^2$。

**边际分布**：$\mathbf y_g\sim\mathcal N(\mathbf X\boldsymbol\beta_g,\mathbf V_g)$。REML 通过最大化

$$
\ell_R(\sigma_g^2,\mathbf D_g^*) =
-\tfrac12\left[\log|\mathbf V_g|+\log|\mathbf X^\top\mathbf V_g^{-1}\mathbf X|
+\mathbf r_g^\top\mathbf V_g^{-1}\mathbf r_g\right]+\text{const}
$$

来估方差成分；给定方差成分后，固定效应的 GLS 估计为

$$
\hat{\boldsymbol\beta}_g=(\mathbf X^\top\mathbf V_g^{-1}\mathbf X)^{-1}\mathbf X^\top\mathbf V_g^{-1}\mathbf y_g,
\qquad
\operatorname{Var}(\hat{\boldsymbol\beta}_g)=(\mathbf X^\top\mathbf V_g^{-1}\mathbf X)^{-1}.
$$

**关键点**：随机效应不是"多几个参数"，而是通过 $\mathbf D_g^*$ 把同类样本"收缩"到一起；方差成分由数据估，参数个数不随 batch 数增长。

### 2.2 为什么不用固定效应建模高维批次

设 batch 有 $B$ 个水平。固定效应模型要估 $B-1$ 个哑变量：

- **自由度损失**：每个基因多消耗 $B-1$ 个 df；当 $B$ 与样本量同阶（如数十个 pool、每 pool 6 个供体），残差 df 所剩无几。
- **效力下降甚至反向**：Pinheiro & Bates（原文引文 31）指出，当批次数相对样本量很大时，把批次放进固定效应可能**比不放更差**（估计的 batch 效应本身噪声大，且与处理变量共线）。
- **随机效应**：假设 batch 效应来自 $\mathcal N(0,\tau^2)$，只估一个 $\tau^2$，用收缩（BLUP）吸收；既不耗 df，又能控制 FPR。原文 Supplementary Fig. 3 正是这一点的实证。

### 2.3 为什么不用负二项 GLMM

- **计算**：每个基因要数值积分/拉普拉斯近似，数百万次拟合不可行（原文 Fig. 2C：GLMM 比 dreamlet 慢一个数量级以上）。
- **有限样本校准**：GLMM 的 Wald/F 检验依赖渐近正态；在 25–75 个受试者时，系数估计的零分布偏离正态，FPR 膨胀（原文 Fig. 2E，MAST/glmer 在小样本时 FPR 明显 >5%）。
- **权重 + 线性模型**：Law 等（voom）已经证明，对 log 变换计数做**加权线性回归**，在有限样本下比 NB 回归的假设检验更稳；dream 把这一结论推广到 LMM。

> ⚠️ **导师批注**：这不是说 LMM 在生物学上"更真"，而是说它在**有限样本的频率学性质**上更可控。代价是 log 变换的方差不稳定（靠权重修正）与计数分布的近似。写讨论时必须说清这是"工程-统计折中"，不是"生成模型的胜利"。

### 2.4 LMM 的残差自由度：$\chi^2$ 混合分布匹配

线性模型里 $d_r=n-p$ 是精确的。LMM 里残差

$$
\mathbf r=(\mathbf I-\mathbf H)\mathbf y,\qquad
\mathbf H=\mathbf X(\mathbf X^\top\mathbf V^{-1}\mathbf X)^{-1}\mathbf X^\top\mathbf V^{-1}
\;(\text{固定效应 hat matrix})
$$

但更完整的 hat matrix 还包含随机效应部分。可以证明（原文 Supplementary Information）：

$$
\frac{\mathbf r^\top\mathbf V^{-1}\mathbf r}{\sigma^2}
\;\sim\;\sum_{i=1}^{n}\lambda_i\chi^2_1,
$$

其中 $\lambda_i$ 是 $(\mathbf I-\mathbf H)\mathbf V$（或其对称化）的特征值，**依赖数据与估计出的方差成分**。原文的做法：

> "We match the expected value of this mixture distribution using a single chi-square and use its degrees of freedom to approximate the $df_r$."

因为 $\mathbb E[\chi^2_d]=d$，匹配均值即 $\hat d_r=\sum_i\lambda_i$。**注意 $\mathbf H$ 在 LMM 下不是幂等矩阵**，所以 $\sum_i\lambda_i\neq n-\operatorname{tr}(\mathbf H)$；若把 $\sum_i\lambda_i$ 记作 $\operatorname{tr}((\mathbf I-\mathbf H)^2)$（$\mathbf H$ 对称化后），则

$$
\boxed{\;\hat d_r=\operatorname{tr}\!\left[(\mathbf I-\mathbf H)^2\right]
=n-2\operatorname{tr}(\mathbf H)+\operatorname{tr}(\mathbf H^2)\;}
$$

- 线性模型：$\mathbf H$ 幂等，$\operatorname{tr}(\mathbf H^2)=\operatorname{tr}(\mathbf H)=p$ → $d_r=n-p$（精确）。
- LMM：$\operatorname{tr}(\mathbf H^2)\le\operatorname{tr}(\mathbf H)$，故 $d_r$ 介于 $n-p$ 与 $n-\operatorname{tr}(\mathbf H)$ 之间；随机效应越强，$\mathbf H$ 越"收缩"，$d_r$ 越大（有效参数更少）。
- **近似性**：用单一 $\chi^2$ 匹配混合分布只匹配了一阶矩；更高阶矩不保证，样本量越大越准。这就是原文说"exact for linear models, approximate for LMM"的含义。

**实现层面**（务必区分三个量）：

| 量 | 来源 | 用途 |
|---|---|---|
| `fit$df.residual` | `lmerTest` 的 **contest 检验自由度** | `summary()` 的 p 值，**不是** EB 用的 df |
| `fit$rdf` | `variancePartition:::rdf.merMod()`：`n − 2tr(H) + tr(H²)` 的低秩实现 | **EB 的 `d_g`**（`eBayes,MArrayLM2` 临时把 `df.residual` 换成 `rdf`） |
| `n − p` | 普通线性模型 | 仅作敏感性对照 |

> ⚠️ **导师批注**：这是本轮审核纠正的最重要一处。手算 moderated t 时若用 `fit$df.residual`，会与 `eBayes` 的结果对不上；必须用 `fit$rdf`。不同实现（Satterthwaite、Kenward–Roger）会给不同的 $d_r$，**本项目固定用 `rdf.merMod` 的线性近似**（D8），并在敏感性分析里与 $n-p$ 对比。

### 2.5 经验贝叶斯 moderated t（Smyth 2004）

**问题**：每个基因独立估 $s_g^2$，在 $n$ 小时极不稳定；方差大的基因 t 小、方差小的基因 t 大，排序被噪声主导。

**先验**：假设

$$
\frac{1}{s_g^2}\sim\frac{1}{d_0 s_0^2}\chi^2_{d_0}
\quad\Longleftrightarrow\quad
s_g^2\sim\text{Inv-}\Gamma\!\left(\frac{d_0}{2},\frac{d_0s_0^2}{2}\right).
$$

**后验**（共轭）：

$$
\tilde s_g^2=\frac{d_0s_0^2+d_gs_g^2}{d_0+d_g},
\qquad
\tilde d_g=d_0+d_g.
$$

**moderated t**：

$$
\tilde t_g=\frac{\hat\beta_g}{\tilde s_g\sqrt{v_g}},
$$

其中 $v_g$ 是 $(\mathbf X^\top\mathbf V^{-1}\mathbf X)^{-1}$ 中对应该系数的对角元。

**零分布（务必注意）**：在线性模型下 $\tilde t_g\sim t_{\tilde d_g}$；但在 **LMM** 下，Supplement 明确说近似自由度只用于 **EB 收缩**，检验的零分布由 **Satterthwaite 或 Kenward–Roger** 近似给出（Eq. 63–64 一带）。因此：

- 不要手写 `2 * pt(-abs(t), df = d0 + d)` 去复现 `dream` 的 p 值；
- 直接用 `eBayes()`/`topTable()` 返回的 `P.Value`；
- 若必须手算，用 `lmerTest::contest()` 或 Satterthwaite/KR 的 df，而不是 $d_0+d_g$。

**$d_0$ 与 $s_0^2$ 的估计**（limma `fitFDist`）：

1. 对 $\log s_g^2$ 有
$$
\log s_g^2 \sim \mathcal N\!\left(\log s_0^2+\psi\!\left(\tfrac{d_g}{2}\right)-\log\!\left(\tfrac{d_g}{2}\right),\ \psi'\!\left(\tfrac{d_g}{2}\right)\right),
$$
其中 $\psi$ 为 digamma、$\psi'$ 为 trigamma。
2. 用样本方差 $v=\operatorname{Var}(\log s_g^2)$ 解
$$
\psi'\!\left(\tfrac{d_0}{2}\right)=v-\overline{\psi'\!\left(\tfrac{d_g}{2}\right)},
$$
用 `trigammaInverse()` 数值求 $d_0$；再令

$$
e_g=\log s_g^2-\psi\!\left(\tfrac{d_g}{2}\right)+\log\!\left(\tfrac{d_g}{2}\right),
\qquad
\hat s_0^2=\exp\!\left(\overline{e_g}\right),
$$

即把每个基因的 $\log s_g^2$ 先校正到"共同 $d_0$"的尺度再取平均（limma 的实现可能含小样本修正；以上是 `fitFDist` 的主干逻辑）。
3. 若 $v$ 太小（基因间方差比抽样噪声还小），$d_0\to\infty$，退化为"所有基因共用方差"。

**LMM 的扩展**：把 $d_g$ 换成 §2.4 的 $\hat d_r$，$s_g^2$ 换成 LMM 的残差方差估计，$v_g$ 换成 GLS 的方差。这一替换在 LMM 下仍成立（近似）；**公式见 Supplementary Information Eq. 63–64**（MOESM1，公开 30 页）。本项目按 Supplement 的 `tr(A)=n−2tr(H)+tr(HHᵀ)` 与 `variancePartition::rdf.merMod()` 落地（见 D8）。

**F 检验（多系数）**：对 $m$ 个系数的联合检验，

$$
\tilde F=\frac{\hat{\boldsymbol\beta}^\top\mathbf V_\beta^{-1}\hat{\boldsymbol\beta}/m}{\tilde s^2}
\sim F_{m,\tilde d},
$$

用于基因集/通路层面的联合效应（`topTable` 的 `F` 列）。

### 2.6 收缩的直觉

- 一个只有 3 个重复、残差方差估计为 0.8 的基因，和一个有 30 个重复、估计为 0.2 的基因：直接比较会不公平。
- EB 把前者向全局均值 $s_0^2$ 拉：$\tilde s^2=(d_0s_0^2+3\times0.8)/(d_0+3)$。若 $d_0=10$、$s_0^2=0.25$，则 $\tilde s^2=(2.5+2.4)/13=0.377$——比 0.8 小得多。
- 结果：**小样本基因不再因为"方差估大了"而丢掉**，同时 FPR 受控。

---

## ③ 原文逐句拆解

| 原文句 | 中文转述 | 对应实现 | 原文是否给公式 | 缺失项决策 |
|---|---|---|---|---|
| "Use of precision-weighted linear mixed models enables accounting for repeated measures study designs, including technical or biological replicates, high-dimensional batch effects due to sample multiplexing, and variation in sequencing depth and cell number." | LMM + 权重 | `dream()` | 给了定性 | — |
| "Using generalized linear mixed models (GLMM), such as a negative binomial mixed model, can be very computationally demanding, suffer from convergence issues in real data, and produce poorly calibrated p-values on finite samples." | GLMM 的三大问题 | 08 号对照 | 定性 | — |
| "In the case of a linear model with n samples and p covariates (including the intercept), the residual degrees of freedom (df_r) is simply n − p. However, the case of a linear mixed model used here is more complicated." | LMM 的 df 不简单 | §2.4 | 未给公式 | D8：迹近似 |
| "we show that the residual variance estimates follow a distribution given by a weighted mixture of n chi-squared random variables, where the weights depend on both the data and the estimated model parameters." | $\sum\lambda_i\chi^2_1$ | §2.4 | 未给 λ 的形式 | 同上 |
| "We match the expected value of this mixture distribution using a single chi-square and use its degrees of freedom to approximate the df_r." | 匹配均值 | $d_r=\sum\lambda_i$ | 未给推导 | 同上 |
| "this method is exact for linear models, approximate for linear mixed models with any number of random effects, and the approximation improves with sample size." | 精确性边界 | — | — | — |
| "Smyth developed an empirical Bayes approach that borrows information across genes to estimate the residual variance." | EB | `eBayes`/`squeezeVar` | 引文 25 | — |
| "A key value in this calculation is the residual degrees of freedom." | df 是关键 | §2.5 | — | — |

---

## ④ 逐步可运行代码

### 4.1 直接调用 `dream()` / `eBayes()`（理解底层）

```r
# code/04_00_dream_basic.R
# 目的：绕开 dreamlet 封装，看清 dream + eBayes 的两步
suppressPackageStartupMessages({
  library(variancePartition); library(BiocParallel); library(limma)
})

proc <- readRDS("processed/proc_psychad.rds")
ex   <- proc[["Micro_PVM"]]         # EList：E（log2 CPM）+ weights
meta <- as.data.frame(colData(proc))  # ★ 样本元数据在 processed 对象的 colData，不在 ex$targets

form <- ~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + Dx
fit  <- dream(ex, form, data = meta, BPPARAM = SnowParam(8))
fit  <- eBayes(fit)

tt <- topTable(fit, coef = "DxAD", number = Inf, sort.by = "P")
head(tt)
```

> `dream()` 返回 `MArrayLM2`；`eBayes()` 用 §2.5 的 `squeezeVar` 做 EB。`dreamlet()` 只是把这两步按细胞类型批量执行。

### 4.2 手动复算 moderated t（验证理解）

```r
# code/04_01_manual_eb.R
# 取一个基因，手算 s2、s2.tilde、t.tilde，与 eBayes 结果对照
g <- "PTPRG"
i <- match(g, rownames(fit$coefficients))

b   <- fit$coefficients[i, "DxAD"]
v   <- fit$stdev.unscaled[i, "DxAD"]^2
s2  <- fit$sigma[i]^2
d   <- fit$rdf[i]          # ★ EB 用的残差 df；不是 fit$df.residual
d0  <- fit$df.prior
s02 <- fit$s2.prior

s2_tilde <- (d0 * s02 + d * s2) / (d0 + d)
t_tilde  <- b / sqrt(s2_tilde * v)
# ⚠️ LMM 下零分布用 Satterthwaite/KR，不是 t_{d0+d}；下面只是线性模型的近似展示
p_tilde  <- 2 * pt(-abs(t_tilde), df = d0 + d)   # 仅作对照，勿用于报告

cat(sprintf("gene=%s  beta=%.4f  s2=%.4f  d=%.1f  d0=%.2f  s02=%.4f\n", g, b, s2, d, d0, s02))
cat(sprintf("manual: s2_tilde=%.4f  t=%.4f  p(approx)=%.3e\n", s2_tilde, t_tilde, p_tilde))
cat(sprintf("limma : t=%.4f  p=%.3e\n", fit$t[i, "DxAD"], fit$p.value[i, "DxAD"]))
stopifnot(abs(t_tilde - fit$t[i, "DxAD"]) < 1e-6)
```

### 4.3 残差自由度：$n-\operatorname{tr}(\mathbf H)$ 的实现

```r
# code/04_02_df_hat.R
# 对同一基因：rdf.merMod 与 dream 的 fit$rdf 对照
# 说明：本脚本接在 4.1/4.2 之后运行（proc/ex/fit/g 已在会话中）；
# 独立运行请先执行 4.1 的 readRDS 与 4.2 的 g <- "PTPRG"。
suppressPackageStartupMessages({ library(lme4); library(variancePartition) })

# 依赖 4.1 的 proc/ex/fit，以及 4.2 的基因 g
df_ex <- as.data.frame(colData(proc))   # 与 4.1 一致
df_ex$y <- ex$E[g, ]
fit_lmer <- lmer(y ~ scale(Age) + Dx + (1|SubID) + (1|poolID) + (1|Sex),
                 data = df_ex, weights = ex$weights[g, ], REML = TRUE)

# variancePartition 的线性近似：n - 2tr(H) + tr(H^2)
df_rdf <- variancePartition:::rdf.merMod(fit_lmer)
cat("n =", nrow(df_ex), " rdf.merMod =", df_rdf, "\n")
cat("dream fit$rdf =", fit$rdf[i], "  fit$df.residual =", fit$df.residual[i], "\n")
stopifnot(abs(df_rdf - fit$rdf[i]) < 1e-6)
```

> 数值不会完全相等：`dream` 的 df 可能用不同的 hat 近似或对权重做了重标定。**允许差异 <1**，若差很多，检查权重是否传对。

### 4.4 单元测试：已知真值的模拟

```r
# code/04_03_sim_unit_test.R
# 目的：验证 moderated t 的 FPR 与功效（不依赖真实数据）
set.seed(20261018)
suppressPackageStartupMessages({ library(limma); library(lme4) })

sim_one <- function(n_subj = 30, n_rep = 2, n_genes = 2000, n_de = 200,
                    sd_subj = 0.5, sd_resid = 0.3, beta_de = 0.6) {
  subj <- factor(rep(seq_len(n_subj), each = n_rep))
  grp  <- rep(c(0, 1), times = n_subj * n_rep / 2)
  b    <- rnorm(n_subj, 0, sd_subj)
  y <- matrix(rnorm(n_genes * length(subj), 0, sd_resid), n_genes)
  y <- y + b[as.integer(subj)]
  y[seq_len(n_de), ] <- y[seq_len(n_de), ] + beta_de * grp
  list(y = y, grp = grp, subj = subj)
}

d <- sim_one()
design <- model.matrix(~ grp)
fit <- lmFit(d$y, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "grp", number = Inf, sort.by = "none")
cat("FPR (null genes):", mean(tt$P.Value[-seq_len(200)] < 0.05), "\n")
cat("Power (DE genes):", mean(tt$P.Value[seq_len(200)] < 0.05), "\n")
```

**预期**：FPR ≈ 0.05（0.035–0.065），Power 随 `beta_de` 上升。若 FPR 远高于 0.05，先查模拟中的批次/受试者结构是否被模型正确指定。

### 4.5 收敛与失败诊断

```r
# code/04_04_convergence.R
# 查看每个细胞类型失败模型的比例与原因
details(fit)                     # dreamletResult 的 details
seeErrors(fit)                   # 若存在
```

> 失败通常来自：某细胞类型某协变量无变异（全为同一水平）、随机效应方差为 0、完全共线。处理：`dropRedundantTerms()` 自动删项；或手动去掉该细胞类型的该协变量。

---

## ⑤ 参数解释与调参空间

| 参数 | 本项目取值 | 可选 | 改变后果 |
|---|---|---|---|
| `formula`（DE） | 见 05 号 | 增删协变量 | 决定"疾病效应"是否被混杂；**不能用与结局相关的变量做随机效应**（会吸收效应） |
| `REML` | `TRUE`（lme4 默认） | `FALSE` | 方差成分估计略不同；小样本推荐 REML |
| `d0` 估计 | limma `squeezeVar` 默认 | 手动设 `d0` | 手动设会改变收缩强度；不建议 |
| `robust`（eBayes） | `FALSE` | `TRUE` | 稳健 EB 能抗离群基因，但会改变零分布；本复现不用 |
| `trend`（eBayes） | `FALSE` | `TRUE` | 若均值–方差趋势未被 voom 权重完全吸收，可加 `trend=TRUE`；会与 03 号重复校正，慎用 |
| `weights` | 03 号输出 | 1 | 全 1 = 不加权，FPR/功效都会变 |
| `BPPARAM` | `SnowParam(8/16)` | — | 只影响速度 |
| df 近似 | `n − tr(H)` | `n − p`、Satterthwaite | 只影响 p 值的轻微偏移；敏感性分析用 |

---

## ⑥ 预期结果与 QC 检查点

| 检查点 | 预期 | 允许波动 | Go/No-Go |
|---|---|---|---|
| `fit$df.prior`（$d_0$） | 数十到数百（基因间共享） | — | ⚠️ 若 $d_0<1$，说明基因间异质性极大，检查是否混入不同细胞类型 |
| `fit$s2.prior`（$s_0^2$） | 与 `median(fit$sigma^2)` 同量级 | ±50% | ⚠️ |
| 残差 df | 接近 $n-p$ 但略大 | — | ✅ |
| 失败模型比例 | <1% | — | ❌ |
| QQ 图（p 值） | 零假设基因贴合对角线，尾部上扬 | — | ❌ 整体偏离 → 模型/权重错 |
| 模拟 FPR | 0.035–0.065 | — | ❌ |
| 模拟 Power（beta=0.6, n=30） | >0.6 | — | ⚠️ |
| `PTPRG` 的 logFC（AD microglia） | 1.59 | ±10% | ❌ |
| `PTPRG` 的 p | 2.94e-28（原文） | 同数量级 | ⚠️ 受样本版本影响 |

---

## ⑦ 常见陷阱与排错表

| 现象 | 根因 | 解决 |
|---|---|---|
| `dream()` 报 "boundary (singular) fit" | 随机效应方差→0（如某协变量在多数细胞类型无变异） | 该细胞类型去掉该随机效应；或 `dropRedundantTerms()` |
| 所有基因 p 值都极小 | 模型把结局变量放进了随机效应（如 `(1|Dx)`），或权重严重失衡 | 检查 formula；DE 模型里结局必须是固定效应 |
| p 值整体偏大（保守） | $d_0$ 过大 / 权重过度收缩 | 检查 `span` 与 `prior.count`；对照 03 号诊断 |
| `eBayes` 后 t 与手算不一致 | 用了 `robust=TRUE` 或 `trend=TRUE` | 显式设 `FALSE` 复算 |
| `hatvalues` 维度不符 | 有 `NA` 或样本被过滤 | `na.omit()`；确保 `weights` 与 `data` 行对齐 |
| 某细胞类型 df 很小 | 样本少 + 随机效应多 | 合并稀有细胞类型；减少随机效应 |
| 结果与 limma 差很多 | limma 没有随机效应 | 预期内；13 号专门量化 |

---

## ⑧ 与原文图表的对应关系

| 原文图 | 本模块产物 | 说明 |
|---|---|---|
| Fig. 1（LMM 与 EB 两步） | `processed/fit_<dataset>.rds` 的 `df.prior`、`s2.prior` | 示意复刻 |
| Fig. 2D（模拟 AUPR） | `code/04_03_sim_unit_test.R` 的扩展 | 08 号做大规摸 |
| Fig. 2E（置换 FPR） | 同上 | 08 号 |
| Methods「Empirical Bayes」 | 本文件 §2.5 | 公式与代码 |
| Supplementary Information（df 推导） | 本文件 §2.4 | Supplement Eq. 63–64（MOESM1） |

> **交付物**：`fit_<dataset>.rds`（05 号生成）、`logs/04_sessionInfo.txt`、`figures/04_qq_*.png`。完成后进入 `05_差异表达分析流程与结果提取.md`。
