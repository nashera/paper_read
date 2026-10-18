# PROJECT_SPEC —— Dreamlet 大规模单细胞差异表达 复现项目规范（L1）

> **上位规范**：`../REPRO_SPEC.md`（通用规范 **v2.1**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（L0 在本项目的实例化；不得违背 L0 红线）
> **本项目**：Hoffman GE, Lee D, Bendl J, Prashant NM, Hong A, Casey C, Alvia M, Shao Z, Argyriou S, Therrien K, Venkatesh S, Voloudakis G, Haroutunian V, Fullard JF, Roussos P. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* **Nature Communications (2026) 17:9597.** DOI: [10.1038/s41467-026-75680-8](https://doi.org/10.1038/s41467-026-75680-8)
> **项目目录**：`20261018_dreamlet_pseudobulk_DE_large_scale/`
> **启动日期**：2026-10-18　**本文件版本**：**v1.0**（对应 L0 v2.1）
>
> **效力**：L0 `../REPRO_SPEC.md` > **本文件 L1** > 14 号文档 §五 QC 门槛（L2）> 各模块文档（L3）。
> **D 表权威位置**：`14_复现质控清单、排错与批判性评述.md` §八；本文件 §8 只维护索引。

---

## 1. 项目定义

### 1.1 一句话目标

以作者公开的 `dreamlet` R 包与 `dreamlet_analysis` 分析代码为基准，**重跑并核验**"pseudobulk + 两阶段精度权重 + 线性混合模型 + 经验贝叶斯 moderated t"这一统计流水线，在三个公开队列（COVID-19 COMBAT、前列腺癌骨转移 Kfoury、阿尔茨海默病 PsychAD）上逐图重建原文证据链，并补齐原文正文与 Supplementary Information（MOESM1，30 页）未写清的算法细节与参数口径。

### 1.2 原文核心结论（复现时逐条验证）

1. `dreamlet` 以 pseudobulk + **精度加权线性混合模型（LMM）** 对每个细胞簇、每个基因拟合，可建模重复测量（技术/生物学重复）、多重样本混合（multiplexing）造成的高维批次效应、以及测序深度与细胞数差异（Fig. 1）。
2. 两级精度权重：① 由 Poisson 计数模型初始化（细胞数、reads 数越多精度越高）；② 用初始权重估计**经验均值–方差趋势**，得到第二轮权重（Fig. 1；Methods「Two-stage weighting」）。在 LMM 上扩展 Smyth 经验贝叶斯 moderated t，残差自由度用"χ² 混合分布匹配单一 χ²"近似（Methods）。
3. 计算性能：1000 donors × 8 cell types × 4.1 M cells 的 pseudobulk 用时 **31 min、峰值内存 26 GB**（Fig. 2A–B）；12 cell types × 326 subjects 的 DE 用时 **45 CPU-min（wall 15 min）**（Fig. 2C）。
4. 统计性能：在 Cromwell 等（muscat）模拟框架下 AUPR 匹配或优于现有方法，效应量估计最准（Fig. 2D）；在真实 snRNA-seq 上置换疾病标签，只有 dreamlet 在**能建模随机效应**的方法中把 FPR 控制在 5%（Fig. 2E）。
5. COVID-19（674K cells / 110 donors）：classical monocyte 中疾病状态解释的中位方差 11.4%，3250/12,472 个基因 >25%；mild COVID vs healthy 有 **1824 个 DEG**（Fig. 3）。
6. 前列腺癌骨转移（Kfoury，9 例患者）：单核细胞簇中疾病状态对 525/1899（27.6%）基因解释 >10% 方差；tumor vs involved bone marrow **157 个 DEG**（Fig. 4）。
7. AD 队列（586 samples / 299 donors / 1.4 M nuclei / 22 细胞类型）：microglia 中受试者解释的中位方差 38.4%（技术重复可重复性良好）；AD 状态仅 101 个基因 >5%，但 **932 个基因受 pool 影响 >5%**（高维批次效应）；microglia 有 **1037 个 DEG**，**PTPRG log2FC = 1.59、p = 2.94e-28**（Fig. 5–6）。
8. 基因集分析（camera/zenith，full spectrum of t-statistics）：microglia 特异上调 p38MAPK cascade（含 **DUSP10**）；神经元亚型上调突触组装/谷氨酸通路（Fig. 6E–F）。

### 1.3 交付物清单

| 类别 | 文件 | 状态 |
|---|---|---|
| 原文 | `2026_Efficient_differential_expression_analysis_of_large-scale_single-cell_transcriptomics_data_using_Dreamlet.pdf` | ✅ 已归档 |
| 项目规范 | **`PROJECT_SPEC.md`（本文件，L1）** | ✅ v1.0 |
| 通用规范 | `../REPRO_SPEC.md`（L0 v2.1） | ✅ |
| 总览与目录 | `00_目录与复现总览.md` | ✅ |
| 模块文档 | `01`–`13` 共 13 篇 | ✅ |
| 质控与批判 | `14_复现质控清单、排错与批判性评述.md`（NN-1） | ✅ |
| 背景速查 | `15_背景知识精讲与方法学原理速查.md`（NN） | ✅ |
| 代码 | `code/*.R`、`code/*.py`、`code/*.sh`（编号与文档一致） | ⏳ S6 实现 |
| 运行日志 | `logs/DECISIONS.md`、`logs/<NN>_sessionInfo.txt`、`logs/00_checksums.md5` | ⏳ S6 生成 |

### 1.4 模块划分与编号表

| 编号 | 文档标题 | 对应原文位置 | 状态 |
|---|---|---|---|
| 00 | 目录与复现总览（`00_目录与复现总览.md`） | Fig. 1；全文 | ✅ |
| 01 | 环境搭建与数据获取（`01_环境搭建与数据获取.md`） | Data/Code availability；Methods「Software versions」 | ✅ |
| 02 | 数据接入与 pseudobulk 聚合（`02_数据接入与pseudobulk聚合.md`） | Fig. 1；Methods「Efficient computation of pseudobulk」 | ✅ |
| 03 | 两阶段精度权重与均值–方差趋势（`03_两阶段精度权重与均值方差趋势.md`） | Fig. 1、3B、4B、5C；Methods「Two-stage weighting」 | ✅ |
| 04 | 精度加权 LMM 与经验贝叶斯 moderated t（`04_精度加权LMM与经验贝叶斯moderated_t.md`） | Methods「Precision-weighted LMM」「Empirical Bayes shrinkage」 | ✅ |
| 05 | 差异表达分析流程与结果提取（`05_差异表达分析流程与结果提取.md`） | Fig. 3E–F、4F、6A–D | ✅ |
| 06 | 方差分解与批次效应量化（`06_方差分解与批次效应量化.md`） | Fig. 3C–D、4C–E、5D–F | ✅ |
| 07 | 基因集分析：camera 与 zenith（`07_基因集分析camera与zenith.md`） | Fig. 3G–H、4G–H、6E–F | ✅ |
| 08 | 计算与统计性能基准复现（`08_计算与统计性能基准复现.md`） | Fig. 2A–E；Fig. S1–S3 | ✅ |
| 09 | 案例 I：COVID-19 重症单核细胞（`09_案例I_COVID19重症单核细胞.md`） | Fig. 3；Fig. S4 | ✅ |
| 10 | 案例 II：前列腺癌骨转移（`10_案例II_前列腺癌骨转移.md`） | Fig. 4；Fig. S5 | ✅ |
| 11 | 案例 III：AD 队列数据生成与上游处理（`11_案例III_AD队列数据生成与上游处理.md`） | Fig. 5A–B；Fig. S6；Methods 湿实验与 QC | ✅ |
| 12 | 案例 IV：AD 差异表达与下游分析（`12_案例IV_AD差异表达与下游分析.md`） | Fig. 5C–H、Fig. 6、Fig. S7–S8 | ✅ |
| 13 | 方法对照：limma/DESeq2/Harmony 与跨队列一致性（`13_方法对照与跨队列一致性.md`） | Fig. S9–S18 | ✅ |
| 14 | 复现质控清单、排错与批判性评述（`14_复现质控清单、排错与批判性评述.md`，NN-1） | 全文 | ✅ |
| 15 | 背景知识精讲与方法学原理速查（`15_背景知识精讲与方法学原理速查.md`，NN） | 全文 | ✅ |

---

## 2. 目录与命名（本项目的具体化）

- 项目 slug：`dreamlet_pseudobulk_DE_large_scale`；目录 `20261018_dreamlet_pseudobulk_DE_large_scale/`。
- 数据源目录：

| 子目录 | 内容 | 来源 |
|---|---|---|
| `raw/combat_covid19/` | COMBAT COVID-19 单细胞 h5ad + 临床变量表 | COMBAT Consortium（Cell 2022）；数据门户/ArrayExpress（见 01 号 §2.3） |
| `raw/kfoury_prostate/` | Kfoury 前列腺癌骨转移 h5ad | GEO **GSE143791** |
| `raw/psychad_ad/` | PsychAD/NPS-AD snRNA-seq h5ad（public release） | Synapse **syn51188606**（PsychAD_public） |
| `raw/benchmark/` | 性能基准用 H5AD、模拟参数、置换标签 | 由本项目 08 号脚本从上述数据派生 |
| `processed/` | `.rds` 对象（pseudobulk、processed assays、dreamlet fit、VP 表） | 本项目脚本 |
| `results/` | `.tsv/.csv` 结果表（topTable、variance partition、zenith、benchmark） | 本项目脚本 |
| `figures/` | `.pdf` + `.png`（300 dpi） | 本项目脚本 |
| `code/` | `NN_MM_<verb>.R/.py/.sh` | 本项目脚本 |
| `logs/` | sessionInfo、运行日志、`DECISIONS.md`、校验和 | 本项目脚本 |

- 命名特例：
  - 细胞类型标签**一律使用原文/作者代码的字符串**（区分大小写）：AD 用 `Micro_PVM`、`EN_L2_3_IT`、`IN_PVALB` 等 22 类；COVID 用 `major_subset` 原值（`cMono`、`ncMono` 等）；Kfoury 用 `cells` 原值。禁止另造缩写。
  - 基因 ID：**AD 与 Kfoury 的 h5ad 用 Ensembl ID**，作图/基因集分析前统一转 SYMBOL（`org.Hs.eg.db`，版本见 §6）；COMBAT 作者代码直接用 `rowData(sce)$feature_name` 作 SYMBOL。转换丢失清单写入 `logs/`。
  - 图号前缀 `Fig<N><L>_`（如 `Fig2A_pseudobulk_memory.pdf`）；补充图用 `FigS<N><L>_`。

---

## 3. 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `min.counts` | `min.count` | 作者 Kfoury/PsychAD 脚本写作 `min.counts`；`dreamlet::processAssays()` 形参为 `min.count`（源码 `R/processAssays.R`） |
| `fitExtractVarPartModel`（variancePartition 旧接口） | `dreamlet::fitVarPart` | 本项目统一用 dreamlet 封装；底层仍调用 variancePartition |
| `dreamletResult` / `dreamletProcessedData` 混写 | 按对象类型分别书写 | `processAssays()` 返回 `dreamletProcessedData`；`dreamlet()` 返回 `dreamletResult` |
| `AD` / `Alzheimer's disease` / `AD case` | `AD`（分组水平），`Control`（对照） | 全项目统一；`Dx` 因子水平 `c("Control","AD")` |
| `pseudobulk` / `pseudo-bulk` / `pseudo bulk` | `pseudobulk` | L0 §3.3 |
| `10X pool` / `poolID` / `sample pool` | `poolID`（变量名），"10X pool"（叙述） | 作者代码列名 `poolID` |

---

## 4. 产物契约表

| 产物路径 | 生产者(唯一) | 格式/结构（列名·尺度·单位·坐标系） | 消费方 | 生成脚本 |
|---|---|---|---|---|
| `raw/combat_covid19/*.h5ad` | 外部（01） | AnnData；`X` 为原始计数（作者代码 `counts(sce)=assay(sce,"X")`）；`obs` 含 `major_subset`、`scRNASeq_sample_ID`、`donor_id`、`Institute`、`GEX_region` | 02, 09 | 下载脚本 `code/01_01_fetch_data.sh` |
| `raw/kfoury_prostate/*.h5ad` | 外部（01） | AnnData；`obs` 含 `cells`（细胞类型）、`ID`（样本）、`subject.id`、`Status`（Benign/Distal/Involved/Tumor）、`subject.status` | 02, 10 | `code/01_01_fetch_data.sh` |
| `raw/psychad_ad/*.h5ad` | 外部（01） | AnnData；`obs` 含 `SubID`、`Channel`、`poolID`、`Dx_AD`、`Sex`、`Age`、`subclass`、`class`、`round_num` | 02, 12 | `code/01_01_fetch_data.sh` |
| `processed/pb_<dataset>.rds` | 02 | **list**：`$pb`（`SingleCellExperiment`，每个细胞类型一个 assay，genes × samples，值为**原始计数之和**，非负整数）、`$cellCounts`（**样本 × 细胞类型**的细胞数矩阵）、`$config`（cluster_id/sample_id 等）；`$pb` 的行名已按 `id_map` 转为 SYMBOL | 03, 05, 06, 07 | `code/02_01_pseudobulk.R` |
| `processed/id_map_<dataset>.tsv` | 02 | 列：`original_id`、`symbol`、`mapped`（TRUE/FALSE）；记录 Ensembl→SYMBOL 的映射与丢失清单；**pb 的唯一行名口径** | 06, 07, 11, 12 | `code/02_01_pseudobulk.R` |
| `processed/proc_<dataset>.rds` | 03 | `dreamletProcessedData`（list）；每 assay 为 `EList`：`E`=log2 CPM（TMM 归一化后，`augmentPriorCount` 自适应先验计数）；`weights`=精度权重（基因 × 样本，非负）；`targets` 含 `lib.size`；`colData(proc)` 含样本元数据 | 05, 06, 07, 13 | `code/03_01_process_assays.R` |
| `processed/weights_poisson_<dataset>.rds` | 03 | **list**（每个细胞类型一个数值矩阵，genes × samples）；第一轮 Poisson 近似精度权重（未平滑）；行名 = SYMBOL，列名 = 样本 ID | 03, 05 | `code/03_01_process_assays.R` |
| `processed/fit_<dataset>.rds` | 05 | `dreamletResult`（list of `MArrayLM2`）；含 `coefficients`、`stdev.unscaled`、`sigma`、`df.residual`（lmerTest 检验 df）、**`rdf`（EB 用的残差 df）**、`t`、`p.value`、`df.prior`、`s2.prior` | 06, 07, 12, 13 | `code/05_01_dreamlet_fit.R` |
| `processed/sce_psychad_qc.rds` | 11 | `SingleCellExperiment`：QC 后、未注释；行名 SYMBOL；`colData` 含 `SubID`/`Channel`/`poolID`/`Dx`/`Sex`/`Age`/`pass_cell` | 11, 12 | `code/11_03_qc.R` |
| `processed/sce_psychad_annot.rds` | 11 | `SingleCellExperiment`：QC + 22 类 `subclass` 注释；行名 SYMBOL | 12, 13 | `code/11_05_annotate.R` |
| `results/topTable_<dataset>_<coef>.tsv` | 05 | 列：`ID`（SYMBOL）、`logFC`（log2 尺度）、`AveExpr`（log2 CPM）、`t`、`P.Value`（双侧）、`adj.P.Val`（BH，**细胞类型内**）、`adj.P.Val.studywide`（**同一 coef** 下跨全部细胞类型合并 BH）、`B`、`assay`（细胞类型）、`coef`；脚本派生列 `SE = logFC/t`（`topTable()` 本身不返回 SE）。**不再单独产出 study-wide 文件** | 09, 10, 12, 13 | `code/05_02_topTable.R` |
| `results/vp_<dataset>.tsv` | 06 | 列：`gene`、`assay`、各方差来源（`Subject`、`poolID`、`Age`、`Sex`、`Dx`、`Residuals` …，**比例 0–1**，由 `calcVarPart` 归一化，行内和≈1）；长表另存 `vp_<dataset>_long.tsv` | 12, 13 | `code/06_01_variance_partition.R` |
| `results/vp_gc_<dataset>.tsv` | 06 | **每个 assay 一行**：`assay`、`n`、`spearman_rho`、`p`、`FDR`（BH，跨 assay） | 12 | `code/06_05_gc_vs_pool.R` |
| `results/zenith_<dataset>_<geneset>.tsv` | 07 | 列：`assay`、`coef`、`Geneset`、`NGenes`、`delta`（`= (G/m2)(集合内平均 t − 全体平均 t)`）、`se`、`p.less`、`p.greater`、`PValue`（双侧）、`FDR`（BH，跨全部 assay×基因集）、`Direction`；可选派生列 `t = delta/se` | 09, 10, 12 | `code/07_01_zenith.R` |
| `results/benchmark_<metric>.csv` | 08 | 列：`method`、`n_donors`、`n_cells`、`celltype`、`metric`、`value`、`run`、`seconds` | 08, 14 | `code/08_0*_*.R` |
| `figures/Fig<N><L>_<slug>.pdf/.png` | 对应模块 | 矢量 PDF + 300 dpi PNG；图注写明 n、检验、p 形式 | 00, 14 | 各模块脚本 |
| `logs/00_checksums.md5` | 01 | `md5  filename`（raw/ 全量） | 14 | `code/00_checksum.sh` |
| `logs/DECISIONS.md` | 运行时 | 表：日期/决策/触发原因/依据/影响范围/决策人 | 14 | 手工 + 脚本追加 |

> 纪律：生产者唯一；格式即契约；改格式/名 = 改本表 + 通知消费方 + 记 `DECISIONS.md`（L0 §4.2）。

---

## 5. 术语表（唯一口径）

| 术语 | 本项目定义 |
|---|---|
| **pseudobulk** | 对同一 `sample_id × cluster_id` 内的细胞**原始 counts 按基因求和**得到的样本级表达向量；不做任何跨样本归一化（归一化在 `processAssays()` 内完成） |
| **样本（sample）** | 一次 10x 建库（AD：`Channel`）或一次患者取材（Kfoury：`ID`；COVID：`scRNASeq_sample_ID`）。**不是**受试者 |
| **受试者（subject/donor）** | 生物学个体（AD：`SubID`；Kfoury：`subject.id`；COVID：`donor_id`）。重复测量 → 随机效应 |
| **表达基因（expressed gene）** | 通过 `edgeR::filterByExpr(min.count, min.prop=0.4, min.total.count=15)` 的基因；各队列 `min.count`：COVID 3、Kfoury 10、AD 10 |
| **细胞数门槛** | 样本进入某细胞类型分析所需最少细胞数 `min.cells`：COVID 3、Kfoury 10、AD 5 |
| **精度权重（precision weight）** | 观测值方差倒数的相对值；dreamlet 中第二级权重 `w = 1/ŝ⁴`（趋势拟合的是 √SD：`ŝ` 为均值–方差趋势对 √(残差标准差) 的预测值，故 `1/ŝ⁴ = 1/SD²`），第一级为 Poisson 近似权重 |
| **moderated t** | 经验贝叶斯收缩后的 t 统计量 `t̃ = β̂ / (s̃·sqrt(v))`，`s̃²=(d0·s0²+d·s²)/(d0+d)`，零分布 `t_{d0+d}` |
| **残差自由度 d** | LMM 下用 `n − tr(H)` 近似（H 为含随机效应的 hat matrix）；见 04 号 §2.4 |
| **方差分数（variance fraction）** | 某变量（随机效应方差或固定效应拟合值方差）占总方差的比例；`variancePartition` 输出，行内和≈1（变量共线时可 >1） |
| **study-wide FDR** | 对**同一 coef** 下全部细胞类型 × 全部基因的 p 值合并做 BH 校正；与 `topTable()` 默认的"细胞类型内 BH"区分 |
| **DEG** | `study-wide FDR < 0.05` 且 `\|log2FC\|` 不限（原文主图未设 logFC 阈值）；若需效应量筛选另注明 |
| **阳性（`<X>`⁺）** | 本项目**不使用**阳性/阴性细胞标签；细胞类型为原文注释标签，不做二值化 |
| **AD / Control** | `Dx` 因子水平：`Control`（无 AD 病理）、`AD`（CERAD possible/probable/definite + Braak 3–6 + CDR 0.5–5） |
| **COVID 疾病状态** | `Source` 六水平：`COVID_HCW_MILD`、`COVID_MILD`、`COVID_SEV`、`COVID_CRIT`、`Sepsis`、`Healthy`（以作者代码 `coef_array` 为准） |

---

## 6. 环境与版本钉死

| 组件 | 版本 | 出处 | 偏离说明 |
|---|---|---|---|
| R | 4.6.x | 我方选择（Bioconductor 3.23 要求 R 4.6） | 原文未报告 R 版本 |
| Bioconductor | 3.23 | 我方选择（dreamlet 1.10.0 所在 release） | 需与 dreamlet 版本匹配 |
| **dreamlet** | **1.10.0**（Bioconductor 3.23 release） | Bioconductor | 原文用 **1.1.24**（开发版）；API 以 1.10.0 为准，差异见 14 号 D3 |
| variancePartition | ≥1.36.1（Bioc 3.23 为 1.42.0） | Bioconductor | 提供 `dream()`、`voomWithDreamWeights()`、`augmentPriorCount()`、`rdf.merMod()` |
| limma | 3.52.1（原文）/ ≥3.62（我方） | 原文 Methods | 主版本差异已登记 |
| edgeR | 3.38.0（原文）/ ≥4.0（我方） | 原文 Methods | `filterByExpr` 行为一致 |
| DESeq2 | 1.36.0（原文）/ ≥1.38（我方） | 原文 Methods | 仅用于 13 号对照 |
| muscat | 1.11.2（原文）/ ≥1.16（我方） | 原文 Methods | 仅用于 08 号基准 |
| MAST | 1.22.0（原文）/ ≥1.26（我方） | 原文 Methods | 仅用于 08 号基准（单细胞级 GLMM） |
| lme4 | ≥1.1-35 | 依赖 | — |
| zellkonverter | ≥1.14 | Bioconductor | 读 H5AD（DelayedArray 后端） |
| SingleCellExperiment / SummarizedExperiment | 最新 release | Bioconductor | — |
| zenith | ≥1.4 | Bioconductor/GitHub | camera 的 LMM 适配 |
| GSEABase / GO.db / org.Hs.eg.db | 最新 release（记录日期） | Bioconductor | 基因集版本写入 `logs/` |
| MSigDB | v2023.2 或更新（记录日期） | MSigDB | Fig. 3H 的 Hallmark 集 |
| Python | 3.11 | 我方选择 | 仅用于 h5ad 检查/转换与作图辅助 |
| STARsolo | 2.7.9 | 原文 Methods | 上游比对 |
| cellSNP | 1.2.0 | 原文 Methods | 基因型去多重 |
| vireo | 0.5.6 | 原文 Methods | 去多重 |
| QTLtools-mbv | 1.3 | 原文 Methods | 基因型一致性 |
| Scrublet | ≥0.2.3 | 原文 Methods | 去双胞 |
| Harmony | ≥1.2 | 原文 Methods | 批次校正（仅 13 号对照） |
| Leidenalg | ≥0.10 | 原文 Methods | 聚类 |

> 每脚本结束写 `logs/<NN>_sessionInfo.txt`（R）或 `logs/<NN>_pip_freeze.txt`（Python），并记录完整执行命令（L0 §5.6）。

---

## 7. 统计与随机性

- **随机种子常量**：`20261018`
  - R：`set.seed(20261018)`；Python：`random.seed(20261018)` / `np.random.seed(20261018)`。
  - 用于：置换检验（08 号 FPR）、模拟数据（08 号）、UMAP/t-SNE（11、12 号）、Scrublet、vireo、Harmony、Leiden。
- 项目特有统计约定：
  1. **模型层级**：DE 主分析 `dreamlet()`（LMM）；VP 用 `fitVarPart()`；基因集用 `zenith_gsa()`（等价 camera）。
  2. **FDR 双口径**：`topTable()` 的细胞类型内 BH 与 study-wide BH 均输出，主图用 study-wide（与原文一致），表格两者并列。
  3. **效应量**：`logFC` 一律 log2 尺度；`SE = logFC / t`；95% CI = `logFC ± 1.96·SE`。
  4. **置换零模型**：08 号 FPR 分析中，在**每个 pool 内**按原始 AD/Control 比例置换受试者标签（保留批次构成），重复 ≥20 次。
  5. **相关性**：GC 含量 vs pool 方差分数用 Spearman（`ρ`、`p`、`n` 同时报告）。
  6. **并行**：`SnowParam(n)` 并行不影响数值结果，但 `logs/` 记录 `n`；BLAS 多线程可能带来 <1e-10 级浮点差异，声明不追求逐字节一致。

---

## 8. 决策登记表（D 表索引）

> **权威全文见 `14_复现质控清单、排错与批判性评述.md` §八**；本表只列编号与一句话裁决。

| # | 决策项 | 涉及模块 | 裁决(主路径) | 备选/敏感性 | 依据 | ★ |
|---|---|---|---|---|---|---|
| D1 | 三个队列的数据版本与获取路径 | 01, 09–12 | COMBAT 用作者分析的 h5ad（原始计数在 `X`）；Kfoury 用 GSE143791 作者 h5ad；AD 用 PsychAD public release（syn51188606） | 用 GEO/ArrayExpress 原始 fastq 从 STARsolo 重跑（成本高） | 原文 Data availability；作者代码 | ★ |
| D2 | 基因 ID 口径 | 02, 06, 07, 11, 12 | h5ad 保留原始 ID；**pseudobulk 聚合后立即按 `id_map` 统一转 SYMBOL**（`processed/id_map_<dataset>.tsv` 落盘，记录丢失数）；此后全部文档只消费 SYMBOL | 全程 Ensembl + Ensembl 版基因集 | 作者代码 `rownames(sce)=rowData(sce)$feature_name`（仅 COVID）；其他队列口径原文未报告 | ★ |
| D3 | dreamlet 版本 | 全部 | 1.10.0（release） | 1.1.24（原文开发版，不可从 Bioconductor 归档稳定获取） | Bioconductor | |
| D4 | `min.cells` / `min.count` | 05, 09–12 | COVID 3/3；Kfoury 10/10；AD 5/10 | 敏感性：±50% 重跑 | 作者公开代码；原文未报告 | ★ |
| D5 | AD 的两个模型公式 | 03, 12 | **估权重/VP**：`~ (1\|SubID) + (1\|poolID) + (1\|Sex) + scale(Age) + (1\|Dx)`；**DE**：`~ (1\|SubID) + (1\|poolID) + (1\|Sex) + scale(Age) + Dx` | 把 Sex 作固定效应；权重公式去掉 `(1\|Dx)` | 作者 `PsychAD_r0_analysis.Rmd`；原文未报告 | ★ |
| D6 | study-wide FDR 定义 | 05, 09–12 | 同一 coef 跨全部细胞类型合并 BH | 细胞类型内 BH | 原文 "study-wide FDR" 未给形式化定义 | ★ |
| D7 | 两阶段权重实现 | 03, 04 | **两条路径都跑**：A=显式 `weightsList`=Poisson 矩阵（复刻 Supplement Eq.16）；B=`weightsList=NULL`（作者公开代码原样，依赖 1.1.24 内部默认；在 1.10.0 上退化为全 1）。以能复现 headline 数字者为主路径，另一条作敏感性 | `priorWeightsAsCounts=TRUE`（数值不同，不等价） | 原文 Methods + Supplementary §2.2 Eq.16；作者代码未显式传 `weightsList`（版本差异） | ★ |
| D8 | 残差自由度近似 | 04 | `rdf = n − 2·tr(H) + tr(H²)`（`variancePartition:::rdf.merMod` 的线性近似；线性模型退化为 `n−p`）；EB 用 `fit$rdf`，**不用** `fit$df.residual` | Satterthwaite / Kenward-Roger | 原文只给“χ² 混合匹配均值”；具体近似取自 `rdf_functions.R`；公式原文未报告 | ★ |
| D9 | 基因集来源 | 07 | COVID Fig. 3H 用 MSigDB **Hallmark TNFA**；Fig. 3G 与 Kfoury/AD 的 GSA 用 **GO**（作者代码：MHC-II `GO0042613`、chaperone `GO0044183`、p38 `GO1900744`、synapse GO）；两套都跑 | 仅 Hallmark；仅 GO | 原文引文 74（MSigDB）；作者代码用 GO | ★ |
| D10 | COVID 样本选择 | 09 | 每位 donor 取最重症一个样本；仅保留 Oxford 中心样本 | 保留全部样本（Source 与 Institute 混杂） | 作者代码注释 | ★ |
| D11 | Kfoury 对比定义 | 10 | `Tumor − Involved`、`Tumor − Distal`、`Involved − Distal` 三个 contrast | 仅 Tumor vs Involved | 作者代码 | ★ |
| D12 | AD 队列规模差异 | 12 | 以 PsychAD public release 实际样本数为准，逐项对照原文 586/299/1.4M | 若 release 不含全部样本，改用 Snapshot | 原文 vs public release 可能不同 | ★ |
| D13 | 性能基准硬件 | 08 | 本机无法复现 512 GB 集群的绝对耗时 → 报**相对倍数**与内存峰值 | 租用同规格节点 | 原文 Methods 明确 512 GB、1 核 | |
| D14 | 置换次数 | 08 | 每种样本量 ≥20 次置换，报告均值与 95% 区间 | 100 次（更稳但更慢） | 原文未报告置换次数 | |
| D15 | 相关性检验 | 06, 12 | Spearman（GC 含量 vs pool 方差分数） | Pearson | 原文 Fig. 5F 写 Spearman | |
| D16 | AD 细胞级 QC 阈值 | 11 | UMI ≥500、基因 ≥250、线粒体 <5%、scDblFinder（替代 Scrublet）；基因表达 0.05% 与供体 ≥50 细胞按原文 | UMI 300/1000；线粒体 10% | 原文 p.9 给了基因 0.05% 与供体 <50 细胞；UMI/基因/线粒体**只写 thresholding、未给数值**（`原文未报告`） | ★ |
| D17 | Harmony 对照成分数 | 13 | 前 5 个成分 | 10/20 | 原文 "top five components" | |
| D18 | 多重混合的层级与 pool 数 | 11, 14 | 原文为 **6 donors/pool**；`poolID` 的具体数量从 public release 的 `poolID` 列实测，不写死 60 | 若 release 的 pool 定义不同，按实际列 | 原文 p.4 "pooling six donors"；pool 数原文未报告 | ★ |
| D19 | STARsolo 的 feature 类型 | 11 | 主路径用 `--soloFeatures Gene`（与作者代码读 `Solo.out/Gene/raw` 一致）；敏感性分析用 `GeneFull` | `GeneFull`（含内含子） | 原文 Methods 只写 STARsolo；作者代码目录为 `Gene` | ★ |
| D20 | 外部对照数据不可得 | 13 | Mathys 2023 汇总表若无法获取 → 登记 gap 表，用方法间一致性替代；**不伪造** | 从原文补充表/Synapse 导出 | 原文引文 51；数据获取未在正文给出 | ★ |

---

## 9. QC 门槛

- **数值层**：原文报告的关键数值与复现值允许 **±10%**（计数类 ±5%，p 值只比对数量级）；逐项见 14 号 §五。
- **判据层**：✅ 通过 → 推进；⚠️ 小幅偏差 → 排查后记录放行；❌ 显著偏差 → **第二套独立实现**（如 limma+voom 或 muscat）复核，仍不符则入 gap 表。
- **结题层**：gap 表逐模块填「文献值 / 我方复现值 / 是否一致 / 差异原因 / 处理方式」。
- 全局 Go/No-Go（摘要）：
  - G1 pseudobulk 细胞总数与样本数（AD 1.4 M / 586 / 299；COVID 674 K / 110；Kfoury 9 患者）误差 ≤5%。
  - G2 每细胞类型表达基因数在原文量级（AD microglia 12k 量级；COVID classical monocyte 12,472）。
  - G3 PTPRG 在 microglia 的 log2FC 方向为正且 |log2FC| ∈ [1.4, 1.8]；p < 1e-20。
  - G4 置换 FPR：dreamlet 在 5% ± 1.5% 内；limma/voom 显著膨胀（>10%）。
  - G5 计算性能：pseudobulk 内存随样本数增长的**斜率**与原文同量级（绝对秒数允许硬件差异）。

---

## 10. 偏离声明（对 L0 的非红线偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| §8.1 头部四要素 | 00/14/15 号文档不强制"模块对应原文"四要素（00 为总览、14 为质控、15 为背景） | L0 §1.3 自身将 00/NN-1/NN 定义为非模块文档 | 仅这三篇 |
| §6.4 下载校验 | 三个大队列原始文件 >100 GB 时，只对**实际下载使用的 h5ad 子集**做 md5，不逐文件校验收割版 | 存储与时间成本 | `logs/00_checksums.md5` 覆盖范围在文件头注明 |

> L0 §12 诚实性红线不可偏离。

---

## 11. 修订记录

| 版本 | 日期 | 修改内容 | 修改人 |
|---|---|---|---|
| v1.0 | 2026-10-18 | 首版：16 篇文档（00–15）结构、产物契约表、术语表、D 表索引（D1–D20）、环境版本表 | 复现项目组 |
