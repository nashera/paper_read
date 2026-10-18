# PROJECT_SPEC —— CRC 髓系靶向治疗单细胞复现项目规范（L1）

> **上位规范**：`../REPRO_SPEC.md`（通用规范 **v2.1**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（L0 在本项目的实例化；不得违背 L0 红线）
> **本项目**：Zhang L, Li Z, Skrzypczynska KM, Fang Q, Zhang W, O'Brien SA, He Y, Wang L, Zhang Q, Kim A, Gao R, Orf J, Wang T, Sawant D, Kang J, Bhatt D, Lu D, Li CM, Rapaport AS, Perez K, Ye Y, Wang S, Hu X, Ren X, Ouyang W, Shen Z, Egen JG, Zhang Z, Yu X. *Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer.* **Cell. 2020;181(2):442–459.e29.** DOI: [10.1016/j.cell.2020.03.048](https://doi.org/10.1016/j.cell.2020.03.048)
> **项目目录**：`20261018_CRC_myeloid_scRNA_CSF1R_CD40/`（仓库根下）
> **启动日期**：2026-10-18　**本文件版本**：**v1.0**（对应 L0 v2.1）
>
> **效力**：L0 `../REPRO_SPEC.md` > **本文件 L1** > `14` 号文档 §六 QC 门槛（L2）> 各模块文档（L3）。
> **D 表权威位置**：`14_复现质控清单、排错与批判性评述.md` §十一；本文件 §8 只维护索引。
> **当前阶段**：**S3 方案完成、S4 独立审核通过（有条件，见 `logs/14_review_record.md`）**。本套文档是复现方案，**尚未执行 S6 数据分析**；凡涉及"我方结果"的单元格一律留待 S6 填写，不得预填。

---

## 1. 项目定义

### 1.1 一句话目标

以原文公开的人类 CRC 单细胞数据（GEO **GSE146771**）、小鼠处理模型数据（ENA **PRJEB34105**、ArrayExpress **E-MTAB-8832**）以及 TCGA COAD/READ 与 GTEx colon bulk 数据为输入，**逐图重建**"人类 CRC 髓系细胞图谱 → C1QC⁺/SPP1⁺ TAM 二分 → 细胞互作网络 → 小鼠 anti-CSF1R / anti-CD40 机制"这一证据链，并补齐原文 Methods 未写清或写死的参数、阈值与算法实现细节。

### 1.2 原文核心结论（复现时逐条验证）

| # | 原文主张 | 关键数值/图 | 复现判据（详见 14 号 §六） |
|---|---|---|---|
| C1 | 联合 Smart-seq2 与 10x 两种平台可同时获得"基因深度"与"细胞覆盖"，识别 38/36 个免疫细胞簇与 13 个髓系亚群 | 43,817（10x）/10,468（Smart-seq2）细胞；Fig. 1、S1–S3 | G1–G2：细胞数、簇数、marker 复现 |
| C2 | 人类 CRC 存在两类功能二分、发育路径不同的 TAM：C1QC⁺（补体/吞噬/抗原提呈）与 SPP1⁺（血管生成/ECM/促转移） | hM12/hM13；Fig. 2、3A–3D；GSVA | G3–G4：marker 与通路方向一致 |
| C3 | FCN1⁺ 单核样前体经不同 RTM 中间态分化为 C1QC⁺ 与 SPP1⁺ TAM；SPP1⁺ 也可由 NLRP3⁺ RTM 转化 | RNA velocity、URD、PAGA；Fig. 2E–2F、S4E–S4F | G5：轨迹方向与关键中间态 |
| C4 | 联合 scRNA-seq 与 TCGA/GTEx bulk 的"相关性细胞互作网络"显示 TAM 与 cDC 是网络核心节点；SPP1⁺ TAM 连接 CAF/内皮，C1QC⁺ TAM 连接 T 细胞 | 24 个组织富集簇 + 3 个基质簇；Fig. 3E、S5、S8A | G6：节点度数、边方向、Hub 身份 |
| C5 | 配体–受体分析预测 CXCL10–CXCR3（C1QC⁺→T）、SDC2–MMP2、SPP1–ITGAV、FN1–ITGA5（SPP1⁺→基质） | 1,000 次置换；Fig. 3F、S5D–S5F | G7：LR 对特异性与置换 p 值 |
| C6 | 小鼠 MC38/Renca 与人类髓系亚群可跨物种对应；mM15_Macro-Vegfa 对应人 SPP1⁺ TAM，mM11–mM14 对应人 C1QC⁺ TAM | 15 个小鼠髓系簇；Fig. 3G–3H、S6 | G8：相似性矩阵与映射关系 |
| C7 | anti-CSF1R 优先耗竭增殖中的 mM12/mM14， sparing mM15（Vegfa⁺、Cd274⁺、Arg1⁺、贴附血管）；对应人类 C1QC⁺ 高增殖、SPP1⁺ 低增殖；C1QC^low SPP1^high 患者预后最差 | Fig. 4；S7A–S7D；TCGA 生存 | G9–G10：簇丰度变化、增殖分数、HR 方向 |
| C8 | anti-CD40 优先激活 Ccl22⁺ cDC1（mM07），上调 CD80/CD86/IL-12；随后扩增 Bhlhe40⁺ Th1-like CD4 T 与 CD8⁺ 记忆 T、减少耗竭 T；与 anti-PD1 协同 | Fig. 5–7；S6–S7；TCR/STARTRAC | G11–G13：cDC1 丰度、T 细胞动力学、克隆指标 |

### 1.3 交付物清单

| 类别 | 文件 | 状态 |
|---|---|---|
| 原文 | `2020_Zhang_L_Cell_Single-Cell_Analyses_Inform_Mechanisms_of_Myeloid-Targeted_Therapies_in_Colon_Cancer.pdf` | ✅ 已归档 |
| 项目规范 | **`PROJECT_SPEC.md`（本文件，L1）** | ✅ v1.0 |
| 通用规范 | `../REPRO_SPEC.md`（L0 v2.1） | ✅ |
| 总览与目录 | `00_目录与复现总览.md`、`README.md` | ✅ |
| 模块文档 | `01`–`13` 共 13 篇 | ✅ S3 方案 |
| 质控与批判 | `14_复现质控清单、排错与批判性评述.md`（NN-1） | ✅ |
| 背景与精读 | `15_文献精读、生物学背景与方法学原理速查.md`（NN） | ✅ |
| 代码 | `code/NN_MM_*.R`、`code/NN_MM_*.py`、`code/NN_MM_*.sh`（编号与文档一致） | ⏳ S6 实现 |
| 运行日志 | `logs/DECISIONS.md`、`logs/<NN>_sessionInfo.txt`、`logs/00_checksums.md5` | ⏳ S6 生成 |

### 1.4 模块划分与编号表

| 编号 | 文档标题 | 对应原文位置 | 状态 |
|---|---|---|---|
| 00 | 目录与复现总览（`00_目录与复现总览.md`） | 全文；Fig. 1–7 | ✅ |
| 01 | 环境搭建与数据获取（`01_环境搭建与数据获取.md`） | Data/Code availability；STAR Methods | ✅ |
| 02 | 人类 CRC 数据接入与统一预处理（`02_人类CRC数据接入与统一预处理.md`） | Fig. S1；Methods「Single-cell RNA-Seq data processing」「Batch effect」「Doublet」 | ✅ |
| 03 | 无监督聚类、细胞注释与恶性细胞 CNV（`03_无监督聚类注释与恶性细胞CNV.md`） | Fig. 1、S2、S3；Methods「Unsupervised clustering」「Similarity」「CNV」 | ✅ |
| 04 | 髓系亚群与 TAM 二分法（`04_髓系亚群与TAM二分法.md`） | Fig. 2A–2D、3A–3D、S4；Methods「Signature genes」「GSVA」「M1/M2」 | ✅ |
| 05 | 单核–巨噬细胞发育轨迹（`05_单核巨噬细胞发育轨迹.md`） | Fig. 2E–2F、S4D–S4G；Methods「Developmental trajectory」 | ✅ |
| 06 | 细胞互作网络 I：scRNA-seq × TCGA/GTEx 相关性网络（`06_细胞互作网络I_scRNA_bulk相关性网络.md`） | Fig. 3E、S5A–S5C、S8A；Methods「Correlative cell-cell interactions」 | ✅ |
| 07 | 细胞互作网络 II：配体–受体与置换检验（`07_细胞互作网络II_配体受体与置换检验.md`） | Fig. 3F、S5D–S5F；Methods「Significant ligand-receptor pairs」 | ✅ |
| 08 | 小鼠肿瘤模型与跨物种髓系整合（`08_小鼠肿瘤模型与跨物种髓系整合.md`） | Fig. 3G–3H、S6；Methods 小鼠部分 | ✅ |
| 09 | anti-CSF1R：TAM 差异耗竭与耐药（`09_antiCSF1R_TAM差异耗竭与耐药.md`） | Fig. 4、S7A–S7D；Methods 小鼠处理 | ✅ |
| 10 | anti-CD40：cDC1 激活与 CD40–CD40LG 轴（`10_antiCD40_cDC1激活与CD40_CD40LG轴.md`） | Fig. 5、S6G–S6I、S7E；Methods | ✅ |
| 11 | anti-CD40：T 细胞动态、TCR 与 STARTRAC（`11_antiCD40_T细胞动态TCR与STARTRAC.md`） | Fig. 6、7、S7F–S7M；Methods「TCR analysis」 | ✅ |
| 12 | TCGA/GTEx bulk 验证与生存分析（`12_TCGA_GTEx_bulk验证与生存分析.md`） | Fig. 4I、5J、7J；S5C；Methods「Survival analysis」 | ✅ |
| 13 | 实验验证与湿实验 SOP（`13_实验验证与湿实验SOP.md`） | Fig. 3D、4F、5G–5I、6D–6G、7E–7I；Methods 湿实验 | ✅ |
| 14 | 复现质控清单、排错与批判性评述（`14_复现质控清单、排错与批判性评述.md`，NN-1） | 全文 | ✅ |
| 15 | 文献精读、生物学背景与方法学原理速查（`15_文献精读、生物学背景与方法学原理速查.md`，NN） | 全文 | ✅ |

---

## 2. 目录与命名（本项目的具体化）

- 项目 slug：`CRC_myeloid_scRNA_CSF1R_CD40`；目录 `20261018_CRC_myeloid_scRNA_CSF1R_CD40/`。
- 数据源目录：

| 子目录 | 内容 | 来源 |
|---|---|---|
| `raw/human_crc_10x/` | 人类 CRC 10x 3' 表达矩阵、细胞元数据、TCR（若有） | GEO **GSE146771** |
| `raw/human_crc_smartseq2/` | 人类 CRC Smart-seq2 表达矩阵与细胞元数据；含 5,126 个既往 T 细胞 | GEO **GSE146771**；Zhang et al., Nature 2018 |
| `raw/mouse_mc38/` | MC38 anti-CD40 / isotype 的 10x 表达与 TCR 矩阵 | ENA **PRJEB34105** / ArrayExpress **E-MTAB-8832** |
| `raw/mouse_renca/` | Renca anti-CSF1R / isotype 的 10x 表达矩阵 | ENA **PRJEB34105** / ArrayExpress **E-MTAB-8832** |
| `raw/tcga_coadread/` | TCGA COAD + READ 表达矩阵、临床与 MSI 注释 | UCSC Xena（https://xena.ucsc.edu/） |
| `raw/gtex_colon/` | GTEx colon 表达矩阵与样本注释 | UCSC Xena |
| `processed/` | `.rds` / `.h5ad` 中间对象（QC 后、注释后、整合后、轨迹、网络） | 本项目脚本 |
| `results/` | `.csv/.tsv` 结果表（簇 marker、GSVA、轨迹、网络边、LR 对、生存） | 本项目脚本 |
| `figures/` | `.pdf` + `.png`（300 dpi） | 本项目脚本 |
| `code/` | `NN_MM_<verb>.R/.py/.sh` | 本项目脚本 |
| `logs/` | sessionInfo、运行日志、`DECISIONS.md`、校验和 | 本项目脚本 |

- 命名特例：
  - **细胞簇标签一律使用原文/作者字符串**，区分大小写与连字符：`hM12_TAM-C1QC`、`hM13_TAM-SPP1`、`mM15_Macro-Vegfa`、`mT09_CD8_Tem-Ccl5`。禁止另造 `C1QC_TAM`、`SPP1_TAM` 等等价写法（别名表见 §3）。
  - 基因名：人类 **HGNC symbol**（大写），小鼠 **MGI symbol**（首字母大写、其余小写）；跨物种比较前用一标一 ortholog 表转换并记录丢失基因（D20）。
  - 图号前缀：主图 `Fig<N><L>_`（如 `Fig4I_km_C1QC_SPP1.pdf`），补充图 `FigS<N><L>_`。
  - 表格前缀：`TableS<N><L>_`（与原文补充表号对应）。

---

## 3. 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `C1QC_TAM` / `C1QC+ TAM` / `TAM-C1QC`（产物文件名） | `hM12_TAM-C1QC`（人类）；`mM11–mM14`（小鼠对应簇） | 原文 Fig. 2–4；产物文件名只允许 `hM12_TAM-C1QC` |
| `SPP1_TAM` / `SPP1+ TAM` / `TAM-SPP1`（产物文件名） | `hM13_TAM-SPP1`（人类）；`mM15_Macro-Vegfa`（小鼠对应簇） | 原文 Fig. 2–4 |
| `10X` / `10x Genomics` / `103`（OCR 误读） | `10x`（平台名）；`10x Genomics`（公司名） | 原文 PDF OCR 将 `10x` 识别为 `103`，正文引用时统一写 `10x` |
| `anti-CSF1R` / `aCSF1R` / `CSF1R blockade` | `anti-CSF1R`（正文）；`aCSF1R` 仅用于原文图注转写 | 原文 Fig. 4、S7 |
| `anti-CD40` / `aCD40` / `CD40 agonist` | `anti-CD40`（正文）；`aCD40` 仅用于原文图注转写 | 原文 Fig. 5–7、S6–S7 |
| `activated cDC1` / `Ccl22+ cDC1` / `mM07_cDC1-Ccl22` | 人类：`activated cDC1`（表型）；小鼠：`mM07_cDC1-Ccl22`（簇名） | 原文 Fig. 5、S6G–S6I |
| `Th1-like` / `Bhlhe40+ Th1` / `BHLHE40+ Th1-like` | 人类：`BHLHE40+ Th1-like`；小鼠：`mT03_CD4_Th-Bhlhe40` | 原文 Fig. 7；大小写按物种 |
| `MSI` / `MSI-H` / `MSS` | `MSI-H`（≥2 个单核苷酸位点不稳定）、`MSS`（无位点不稳定） | 原文 STAR Methods「Microsatellite instability testing」 |
| `Ro/e` / `RoE` / `ROE` | `Ro/e`（observed/expected ratio） | 原文 Fig. 2C、S3C；Methods「Tissue distribution」 |

> 规则（L0 §3.3）：改名必须**同时**完成 ① 改所有引用处 ② 更新 §4 契约表 ③ 记 `logs/DECISIONS.md`。

---

## 4. 产物契约表

| 产物路径 | 生产者(唯一) | 格式/结构（列名·尺度·单位·坐标系） | 消费方 | 生成脚本 |
|---|---|---|---|---|
| `raw/**` | 外部（01） | 只读；原始矩阵/元数据/临床表；下载后写 `logs/00_checksums.md5` | 02, 08, 12 | `code/01_01_fetch_data.sh` |
| `processed/human_10x_qc.rds` | 02 | `SingleCellExperiment`；`assays$counts` 为原始 UMI（整数）；`colData` 含 `cell_id`、`patient_id`、`tissue`（Tumor/Normal/Blood）、`platform`、`nGene`、`nUMI`、`percent_mito`、`pass_qc`；基因为 HGNC symbol | 03, 04, 06, 07 | `code/02_01_qc_human_10x.R` |
| `processed/human_smartseq2_qc.rds` | 02 | `SingleCellExperiment`；`assays$counts` 为 gene-level counts；`colData` 同 10x，另含 `prior_T_cell`（TRUE/FALSE）；TPM 存 `assays$tpm` | 03, 04, 06, 07 | `code/02_02_qc_smartseq2.R` |
| `processed/human_10x_norm.rds` | 02 | `SingleCellExperiment`；`assays$logtpm` = log2(TPM+1)，已做 pagoda2 batch-factor 校正；`colData$batch_factor` 记录每个细胞的 BF | 03, 04, 05, 06, 07 | `code/02_03_batch_norm.R` |
| `processed/human_smartseq2_norm.rds` | 02 | 同上；`assays$logtpm` = log2(TPM+1) | 03, 04, 05, 06, 07 | `code/02_03_batch_norm.R` |
| `processed/human_annot.rds` | 03 | **10x** `SingleCellExperiment`；`colData$cluster_id`（原文 `hM01…hT18` 字符串）、`major_type`（Myeloid/ILC/T/B/Epithelial/Endothelial/Fibroblast/Malignant）、`tissue`、`patient_id`、`platform="10x"`；恶性判定见 `processed/human_cnv.rds` 的 `colData$is_malignant` | 04, 12 | `code/03_02b_finalize_annot.R` |
| `processed/human_smartseq2_annot.rds` | 03 | **Smart-seq2** `SingleCellExperiment`；字段同上，`platform="Smart-seq2"`，含 `prior_T_cell`；行名与 `human_smartseq2_norm.rds` 完全一致 | 04, 06, 07 | `code/03_02b_finalize_annot.R` |
| `processed/human_cnv.rds` | 03 | `SingleCellExperiment`；`colData$cnv_score`、`cnv_correlation`、`is_malignant`；`assays$cnv` 为 CNVf 矩阵（基因 × 细胞，已按染色体排序） | 03, 04 | `code/03_03_cnv.R` |
| `processed/myeloid_human.rds` | 04 | **10x** `SingleCellExperiment`；仅 `hM01–hM13`；`assays$logtpm`；`colData$cluster_id`、`tissue`、`patient_id`、`platform="10x"` | 04, 05 | `code/04_01_subset_myeloid.R` |
| `processed/myeloid_human_smartseq2.rds` | 04 | **Smart-seq2** `SingleCellExperiment`；仅 `hM01–hM13`；`assays$logtpm`；`colData$cluster_id`、`tissue`、`patient_id`、`platform="Smart-seq2"` | 04 | `code/04_01_subset_myeloid.R` |
| `results/mouse_cd45_counts.tsv` | 08 | 列：`mouse_id`、`study`、`model`、`treatment`、`timepoint`、`tissue`、`n_cd45`（全 CD45⁺ 细胞数） | 09, 10, 11 | `code/08_01_mouse_qc_annot.R` |
| `results/mouse_cd8_counts.tsv` | 08 | 列：`mouse_id`、`study`、`model`、`treatment`、`timepoint`、`tissue`、`n_cd8`（CD8⁺ T 细胞数） | 11 | `code/08_02b_tcell_annot.R` |
| `results/myeloid_markers.tsv` | 04 | 列：`gene`、`cluster_id`、`avg_log2FC`（log2 尺度）、`pct.1`、`pct.2`、`p_val`、`p_val_adj`（BH，簇内）、`FDR_anova`、`Tukey_adj_p`；一行一个 gene×cluster | 04, 14, 15 | `code/04_02_markers.R` |
| `results/tam_gsva.tsv` | 04 | 列：`pathway`、`cluster_id`（hM12/hM13）、`gsva_score`（每细胞）、`t`、`p`、`FDR`（limma moderated t，BH） | 04, 09, 15 | `code/04_04_gsva.R` |
| `processed/trajectory_velocity.rds` | 05 | `list`：`$dm`（diffusion map，细胞 × 2）、`$velocity`（基因 × 细胞）、`$embedding`（细胞 × 2）；细胞 ID 与 `myeloid_human.rds` 一致 | 05, 15 | `code/05_01_velocity.R` |
| `results/paga_edges.tsv` | 05 | 列：`from`、`to`、`weight`（PAGA connectivity，0–1）、`cluster_from`、`cluster_to` | 05, 15 | `code/05_03_paga.py` |
| `results/celltype_signatures.tsv` | 06 | 列：`gene`、`cluster_id`、`p_anova`、`FDR_anova`、`max_abs_log2FC`、`min_tukey_p`、`auc`、`typical_marker`、`is_signature`（TRUE/FALSE）；满足 FDR<0.05、Tukey |log2FC|>1 & adj p<0.01、AUC>0.65、典型 marker 四条 | 06, 07, 12 | `code/06_01_signatures.R` |
| `results/abundance_tcga.tsv` | 06 | 行：TCGA 样本；列：`cluster_id`；值：signature 基因 z-score 均值（无量纲）；含 `cohort`（COAD/READ） | 06, 12 | `code/06_02_abundance.R` |
| `results/abundance_gtex.tsv` | 06 | 行：GTEx colon 样本；列：`cluster_id`；值同上 | 06, 12 | `code/06_02_abundance.R` |
| `results/correlation_network_edges.tsv` | 06 | 列：`source_cluster`、`target_cluster`、`n_top_genes`（13）、`enrichment_z`、`edge`（TRUE/FALSE，z>1.96）、`tissue`（Tumor/Normal） | 06, 15 | `code/06_04_network.R` |
| `results/lr_pairs_significant.tsv` | 07 | 列：`ligand`、`receptor`、`celltype_ligand`、`celltype_receptor`、`interaction_score`（E_L×E_R）、`p_perm`（1,000 次）、`p_adj`（BH）、`specificity` | 07, 15 | `code/07_02_lr_permutation.R` |
| `processed/mouse_annot.rds` | 08 | **Seurat 5 对象**（合并、粗分）；`meta.data` 含 `major_type`（T/NK/B/Myeloid/DC/Stromal/Tumor）、`model`（MC38/Renca）、`treatment`（Isotype/anti-CD40/anti-CSF1R）、`timepoint`（day2/day10/NA）、`tissue`（Tumor/tdLN）、`mouse_id`、`study`；RNA assay 已 `JoinLayers()` | 08 | `code/08_01_mouse_qc_annot.R` |
| `results/human_mouse_similarity.tsv` | 08 | 列：`human_cluster`、`mouse_cluster`、`mean_logit`、`probability`（0–1）、`model`（MC38/Renca） | 08, 15 | `code/08_03_similarity.R` |
| `results/csf1r_abundance.tsv` | 09 | 列：`cluster_id`、`treatment`、`mouse_id`、`n_cells`、`pct_of_CD45`（0–100）、`study`（Study1/Study2） | 09, 14 | `code/09_01_abundance.R` |
| `results/csf1r_de.tsv` | 09 | 列：`gene`、`cluster_id`、`log2FC`（anti-CSF1R vs isotype）、`p`、`FDR`、`direction` | 09, 15 | `code/09_02_de.R` |
| `results/cd40_dc_abundance.tsv` | 10 | 列：`cluster_id`、`treatment`、`timepoint`、`mouse_id`、`pct_of_CD45` | 10, 14 | `code/10_01_dc_abundance.R` |
| `results/cd40_t_abundance.tsv` | 11 | 列：`cluster_id`、`treatment`、`timepoint`、`tissue`、`mouse_id`、`pct_of_CD8` 或 `pct_of_CD45` | 11, 14 | `code/11_01_t_abundance.R` |
| `results/startrac_indices.tsv` | 11 | 列：`mouse_id`、`treatment`、`timepoint`、`cluster_id`、`expansion`、`migration`、`n_cells`、`n_clones`（无量纲）；transition 另见 `results/startrac_pairwise_transition.tsv` | 11, 15 | `code/11_03_startrac.R` |
| `results/survival_cox.tsv` | 12 | 列：`signature`、`group`、`n`、`HR`、`HR_low`、`HR_high`、`p`、`covariates`（stage+MSI）、`cohort` | 12, 14 | `code/12_02_survival.R` |
| `figures/Fig<N><L>_<slug>.pdf/.png` | 对应模块 | 矢量 PDF + 300 dpi PNG；图注写明 n、检验、p 形式 | 00, 14 | 各模块脚本 |
| `logs/00_checksums.md5` | 01 | `md5  filename`（`raw/` 全量） | 14 | `code/00_checksum.sh` |
| `logs/DECISIONS.md` | 运行时 | 表：日期/决策/触发原因/依据/影响范围/决策人 | 14 | 手工 + 脚本追加 |

**跨模块共享产物补充表（S3 文档落地后追加；与上表同效力）**

| 产物路径 | 生产者(唯一) | 格式/结构 | 消费方 | 生成脚本 |
|---|---|---|---|---|
| `results/ro_e.tsv` | 03 | 列：`cluster_id`、`tissue`、`observed`、`expected`、`ro_e`（无量纲，>1 富集）、`platform`（10x/Smart-seq2） | 06 | `code/03_04_ro_e.R` |
| `results/correlation_gene_cluster.tsv` | 06 | 长表：`gene`、`cluster_id`、`r`（Pearson）、`p`、`n`、`is_self`（TRUE/FALSE） | 06, 12 | `code/06_03_correlation.R` |
| `results/activated_dc_signature_genes.tsv` | 10 | 列：`gene`、`source`（`TableS4` 或 `provisional_8_genes`） | 12 | `code/10_04_export_signature.R` |
| `results/activated_dc_signature_used.tsv` | 12 | 列：`gene`（最终用于打分的基因，已与 bulk 交集） | 12, 14 | `code/12_01_bulk_signature.R` |
| `results/th1_like_signature.tsv` | 12 | 列：`gene`（Th1-like signature，文献/Table S4） | 12, 14 | `code/12_01_bulk_signature.R` |
| `results/bulk_signature_scores.tsv` | 12 | 宽表：`sample`、`C1QC`、`SPP1`、`activated_DC`、`Th1_like`（每列 = signature 基因的 z-score 均值）；GTEx 版为 `results/bulk_signature_scores_gtex.tsv`，列相同 | 12 | `code/12_01_bulk_signature.R` |
| `results/survival_groups.tsv` | 12 | 列：`C1QC_SPP1`（四组标签）、`n`（每组样本数） | 12, 14 | `code/12_02_survival.R` |
| `results/th1_cdc1_correlation.tsv` | 12 | 列：`celltype`（BATF3/activated cDC1/cDC2）、`r`、`p`、`n` | 12, 14 | `code/12_03_correlation_msi.R` |
| `results/cdc1_msi.tsv` | 12 | 列：`comparison`、`n_MSI_H`、`n_MSS`、`median_MSI_H`、`median_MSS`、`p`（Wilcoxon） | 12, 14 | `code/12_03_correlation_msi.R` |
| `results/tam_deg.tsv` | 04 | 列：`gene`、`log2FC`（hM12 vs hM13）、`p`、`FDR`、`direction` | 04, 06, 15 | `code/04_03_deg_tam.R` |
| `results/m1_m2_proliferation.tsv` | 04 | 列：`cluster_id`、`M1_score`、`M2_score`、`proliferation_score`（均为基因集均值） | 04, 09, 15 | `code/04_05_proliferation_m1m2.R` |
| `results/proliferation_scores.tsv` | 09 | 列：`cell_id`、`cluster_id`、`score`（23-gene mean）、`species` | 09, 12 | `code/09_03_proliferation.R` |
| `results/mouse_myeloid_markers.tsv` | 08 | 列：`gene`、`cluster_id`（mM01–mM15）、`avg_log2FC`、`p_val_adj` | 08–11 | `code/08_01_mouse_qc_annot.R` |
| `results/tcr_clonotypes.tsv` | 11 | 列：`cell_id`、`sample_id`、`barcode`、`clone_id`、`cdr3.TRA`、`cdr3.TRB`、`umis.TRA`、`umis.TRB`、`cluster_id`、`tissue`、`treatment`、`timepoint`、`mouse_id` | 11 | `code/11_02_tcr_process.R` |
| `results/startrac_pairwise_transition.tsv` | 11 | 列：`from`、`to`、`transition_index`、`p`、`n_mice` | 11, 15 | `code/11_03_startrac.R` |
| `results/validation_flow_summary.tsv` | 13 | 列：`gate`、`group`、`n`、`mean`、`sd`、`p` | 13, 14 | `code/13_01_flow_gating.R` |
| `results/module_metrics.tsv` | 02–13 | 列：`metric`、`observed`、`module`（长表；各模块追加） | 14 | 各模块脚本 `append_metric()` |
| `results/qc_summary.tsv` | 14 | 列：`metric`、`literature`、`observed`、`abs_diff`、`pct_diff`、`status` | 14 | `code/14_01_qc_summary.R` |
| `results/gap_table.tsv` | 14 | 列：`module`、`metric`、`literature`、`observed`、`consistent`、`reason`、`action` | 14 | `code/14_01_qc_summary.R` |
| `processed/signature_inputs.rds` | 06 | list：`expr`（Smart-seq2 logTPM）、`meta`、`signatures` | 06 | `code/06_01_signatures.R` |
| `processed/correlation_matrices.rds` | 06 | list：`r_gene_cluster`、`abundance_tcga`、`abundance_gtex` | 06, 12 | `code/06_03_correlation.R` |
| `processed/mouse_myeloid_integrated.rds` | 08 | **SingleCellExperiment**；整合后的 15 个髓系簇；`assays$logcounts`；`colData` 含 `cluster_id`（`mM01…mM15`）、`model`、`treatment`、`timepoint`、`tissue`、`mouse_id`、`study` | 09, 10 | `code/08_02_integrate.R` |
| `processed/mouse_t_annot.rds` | 08 | **SingleCellExperiment**；T 细胞子集；`assays$logcounts`；`colData` 含 `cluster_id`（`mT01…mT12`）、`model`、`treatment`、`timepoint`、`tissue`、`mouse_id`、`study` | 11 | `code/08_02b_tcell_annot.R` |
| `processed/myeloid_human.h5ad` | 05 | AnnData；`X` 为 spliced counts，`layers["unspliced"]`；用于 scVelo | 05 | `code/05_01_velocity.R` |
| `processed/urd_object.rds` | 05 | URD 对象；含 pseudotime、tree、segment | 05 | `code/05_02_urd.R` |

> 纪律：生产者唯一；格式即契约；改格式/名 = 改本表 + 通知消费方 + 记 `DECISIONS.md`（L0 §4.2）。

---

## 5. 术语表（唯一口径）

| 术语 | 本项目定义 |
|---|---|
| **TAM（tumor-associated macrophage）** | 原文中以 `hM12_TAM-C1QC`、`hM13_TAM-SPP1` 表示的两类肿瘤富集巨噬细胞；小鼠对应 `mM11–mM15` 中的巨噬细胞簇。**不把 RTM 混称 TAM** |
| **RTM（resident tissue macrophage）** | 正常黏膜富集或相当的巨噬细胞：`hM08_Macro-NLRP3`、`hM09_Macro-PLTP`、`hM10_Macro-IL1B` |
| **C1QC⁺ TAM** | `hM12_TAM-C1QC`；marker 规则：`C1QA/B/C` 平均表达 > 1（log2 TPM）且 >20% 细胞表达，`SPP1` 不满足该规则 |
| **SPP1⁺ TAM** | `hM13_TAM-SPP1`；marker 规则：`SPP1`、`MARCO`、`VEGFA` 平均表达 > 1 且 >20% 细胞表达，`C1QC` 不满足该规则 |
| **组织富集（tissue-enriched）** | `Ro/e > 1`（observed/expected，期望由卡方检验给出）；`Ro/e > 2` 记为强富集（原文 Fig. S8B 用 `+++`） |
| **signature gene（细胞亚群特异基因）** | 同时满足：① ANOVA FDR < 0.05；② 任一簇对 Tukey HSD |log2FC| > 1 且 adj p < 0.01；③ AUC > 0.65；④ 文献典型 marker。四条缺一不可 |
| **self-expressed gene** | 在 Smart-seq2 中该细胞类型平均表达 > 1（log2 TPM）且表达细胞比例 > 20%；在相关性矩阵中将其相关系数置 0 |
| **expressed ligand/receptor** | 某细胞类型中 >20% 细胞满足 log2(TPM+1) > 3；不满足者表达值置 0 |
| **LR interaction score** | `Int(A/B) = E_A(ligand) × E_B(receptor)`；A 为配体细胞类型，B 为受体细胞类型 |
| **M1/M2 score** | 原文 14 个 M1 基因、27 个 M2 基因的平均表达（基因列表见 04 号 §2.6）；仅作对照，不作为主分型 |
| **proliferation score** | 23 个增殖基因（`AURKA…ZWINT`）的平均表达；基因列表见 04 号 §2.6 |
| **MSI-H** | 7 个位点（5 单核苷酸 + 2 五核苷酸）中 **≥2 个单核苷酸位点**不稳定；MSS = 无位点不稳定 |
| **高/低表达分组（survival）** | 按 signature mean expression 的 **55th/45th 分位数**分组：≥55th = high，≤45th = low，中间 10% 不纳入主分析；敏感性分析用中位数二分 |
| **activated cDC1** | 人类：`LAMP3⁺CCR7⁺` 表型；小鼠：`mM07_cDC1-Ccl22`（`Ccl22` 高、`Clec9a/Xcr1` 低）；判定以 cluster 标签 + marker 规则为准 |
| **Th1-like cell** | 人类：`BHLHE40⁺` CD4 T，`IFNG` 高；小鼠：`mT03_CD4_Th-Bhlhe40`，`Bhlhe40` 高、`Ifng` 高 |
| **阳性（`<X>`⁺）** | 流式：按原文图注的门控策略；scRNA：若未特别说明，指该基因在细胞中 UMI > 0；cluster 水平指平均表达 > 1 且表达比例 > 20% |
| **differential depletion（差异耗竭）** | anti-CSF1R vs isotype 的簇丰度（占 CD45⁺ 百分比）变化，单侧 t 检验 p < 0.05；原文未给效应量阈值，另报 fold change |
| **clonal expansion（高克隆）** | 同一 dominant TCR αβ 对在 `Ccl5⁺ Tem` 簇内 ≥3 个细胞为 high clonal，其余为 low clonal（原文 Fig. 6H） |

> 任何「阳性 (<X>⁺)」标签必须给出判定规则；上表为全项目唯一口径（L0 §8.3-5）。

---

## 6. 环境与版本钉死

| 组件 | 版本 | 出处 | 偏离说明 |
|---|---|---|---|
| R | 4.5.1 | 我方选择（Bioconductor 3.22） | 原文未报告 R 版本 |
| Bioconductor | 3.22 | 我方选择 | 需与 Seurat/GSVA 版本匹配 |
| **Seurat** | **5.1.0** | 原文用 2.3.0（人类）/3.1.0（小鼠整合）；我方用 5.1.0 | 大版本差异已登记 D2；`FindClusters`/`IntegrateData` API 已变化 |
| SeuratObject | 5.0.2 | CRAN | — |
| glmnet | 4.1-8 | CRAN；原文 2.0.16 | `cv.glmnet` 接口兼容 |
| GSVA | 1.54.0 | Bioconductor；原文 1.3.0 | `gsva()` 默认参数变化，需显式设 `kcdf="Gaussian"`/`mx.diff=TRUE`（按原文） |
| limma | 3.62.0 | Bioconductor；原文 3.37.11 | — |
| survival | 3.8-3 | CRAN | — |
| survminer | 0.4.9 | CRAN | `ggsurvplot` |
| ComplexHeatmap | 2.22.0 | Bioconductor | — |
| pagoda2 | 1.0.12（GitHub，commit 记录于 `logs/`） | 原文用 pagoda2 batch factor | 若安装失败，主路径改用项目内 R 实现（`code/02_03_batch_norm.R` 附公式），敏感性用 `sva::ComBat` |
| Scrublet | 0.2.3 | Python；原文引用 | — |
| scanpy | 1.10.3 | Python；原文 1.4.4 | PAGA/UMAP API 变化 |
| anndata | 0.10.9 | Python | — |
| scVelo | 0.3.2 | Python | 替代原文 velocyto R 嵌入；速度矢量仍由 velocyto.py 生成 |
| velocyto.py | 0.17.17 | Python；原文 0.17 | — |
| destiny | 3.18.0 | Bioconductor；原文 2.11.3 | diffusion map |
| URD | 1.0.1 | GitHub；原文 1.0.1 | 仅 R 4.2 容器可运行；主路径用 PAGA + scVelo，URD 作敏感性（D10） |
| STARTRAC | GitHub `Japrin/STARTRAC`（commit 记录） | 原文 | R 脚本 |
| Cell Ranger | 2.1.0（原文）；7.2.0（可选重跑） | 10x Genomics | 若仅用 processed matrix 则不执行 |
| kallisto | 0.43.1（原文）；0.50.1（可选重跑） | 原文 | Smart-seq2 原始 FASTQ 不可得（人类隐私），仅用作者 processed counts |
| tximport | 1.9.12（原文）；1.34.0（我方） | Bioconductor | — |
| infercnv | 1.20.0（可选） | Bioconductor | 仅作 CNV 敏感性分析；主路径用原文公式 |
| Python | 3.11.9 | 我方选择 | — |
| FANTOM5 LR pairs | Ramilowski et al. 2015（下载日期记录于 `logs/`） | 原文 | — |
| CellPhoneDB LR pairs | v2.0.0（原文引用 2018 版） | 原文 | 记录版本与下载日期 |
| MSigDB | v2024.1.H（或执行日最新版） | 原文未给版本 | GSVA 主路径用 MSigDB C2 CP + Hallmark；GO BP 作敏感性（D11） |
| GO | Bioconductor `org.Hs.eg.db` 3.19 / `GO.db` 3.19 | 原文用 GSEABase | 记录版本 |
| TCGA/GTEx | UCSC Xena 下载日期 2026-10-18 | 原文 | 记录 `TOIL` 版本与样本过滤 |

> 每脚本结束写 `logs/<NN>_sessionInfo.txt`（R）或 `logs/<NN>_pip_freeze.txt`（Python），并记录完整执行命令（L0 §5.6）。
>
> ⚠️ **本机实测版本说明（2026-10-18）**：当前工作站为 **R 4.4.1 / Python 3.12.4**；上表是按目标环境钉死的 **R 4.5.1 / Bioc 3.22 / Python 3.11.9**。S6 执行前必须二选一：① 用 `renv`/conda 安装目标版本；② 整体降到 R 4.4.1 + Bioc 3.19 + Python 3.12，并在 `logs/DECISIONS.md` 记录。**不得混用**（Seurat/GSVA/Bioc 的年度 release 不兼容）。

---

## 7. 统计与随机性

- **随机种子常量**：`20261018`
  - R：`set.seed(20261018)`；Python：`random.seed(20261018)` / `np.random.seed(20261018)`。
  - 用于：t-SNE/UMAP、Seurat 聚类、Scrublet、RNA velocity 随机采样、PAGA、1,000 次 LR 置换、logistic regression CV、bootstrap。
- 项目特有统计约定：
  1. **差异表达**：原文主路径为 ANOVA + Tukey HSD（所有簇两两比较）；TAM 二分为 limma moderated t。两者均报 `p` 与 BH `FDR`；原文只写 `P < 0.05` 处，双口径（raw/adj）并列。
  2. **GSVA**：按细胞计算 pathway score，再用 limma moderated t 比较 hM12 vs hM13；报 `t`、`p`、`BH FDR`。
  3. **相关性网络**：Pearson（原文明确用 Pearson）；必须同时报 `r`、`p`、`n`（TCGA 样本数）。自表达基因置 0 的规则固定。
  4. **LR 置换**：1,000 次打乱 cluster label；p 值用置换分布的正态近似（原文），另报经验 p（`(b+1)/(B+1)`）作敏感性。
  5. **生存分析**：Cox 多因素校正 `stage` 与 `MSI`；HR 必附 95% CI 与 n；KM 用 `ggsurvplot`；主分析 55/45 分位，敏感性用中位数二分。禁止用最优截断点。
  6. **小鼠丰度**：以 `pct_of_CD45` 为单位；组间比较用单侧 t 检验（原文），同时报双侧 p 与 fold change；n = 生物学重复（小鼠），不是细胞数。
  7. **多重比较**：全项目默认 BH(FDR)；置换检验的 `p_adj` 在 LR 对层面做 BH。
  8. **不可固定项**：Seurat/UMAP 在不同 BLAS/并行下可能有 <1e-6 的嵌入差异；报告中声明，不追求逐字节一致。

---

## 8. 决策登记表（D 表索引）

> **权威全文见 `14_复现质控清单、排错与批判性评述.md` §十一**；本表只列编号与一句话裁决。

| # | 决策项 | 涉及模块 | 裁决(主路径) | 备选/敏感性 | 依据 | ★ |
|---|---|---|---|---|---|---|
| D1 | 数据版本与获取路径 | 01, 02, 08 | 人类用 GSE146771 processed matrix；小鼠用 E-MTAB-8832/PRJEB34105 processed matrix；TCGA/GTEx 用 UCSC Xena | 从 FASTQ 重跑（人类 FASTQ 不可得） | 原文 Data availability；隐私限制 | ★ |
| D2 | 软件大版本 | 全部 | Seurat 5.1.0 + scanpy 1.10.3；原文版本写入 15 号对照表 | 用 renv 锁 Seurat 3.1.0（与原文更近，但安装困难） | 原文 Methods；L0 §5.4 | ★ |
| D3 | 基因组与基因 ID | 02, 03, 08 | 人类 GRCh38/hg38 + HGNC symbol；小鼠 GRCm38/mm10 + MGI symbol；线粒体比例阈值 10%（CD45⁺）/70%（CD45⁻） | GRCh37；Ensembl ID | 原文 Methods | ★ |
| D4 | 细胞 QC 阈值 | 02 | Smart-seq2：library size 与基因数 < median − 3×MAD 剔除；10x：600–4,000 genes、1,600–25,000 UMIs；小鼠沿用 10x 阈值（原文未给） | 放宽/收紧 ±20% 作敏感性 | 原文 Methods；小鼠阈值原文未报告 | ★ |
| D5 | 批次校正 | 02 | pagoda2 batch factor（原文公式）；Smart-seq2 与 10x 分别校正 | `sva::ComBat`；按患者分组建模 | 原文 Methods「Batch effect correction」 | ★ |
| D6 | 双细胞去除 | 02 | Scrublet 默认阈值 + 删除 `CD14⁺IL32⁺` 双细胞簇 | scDblFinder；DoubletFinder | 原文 Methods | ★ |
| D7 | 聚类参数 | 03 | Seurat 5：2,000 HVG、10–25 PCs、resolution 0.6–2.0；逐数据集独立聚类 | 固定 resolution 1.0；Leiden | 原文 Methods | ★ |
| D8 | 细胞注释 | 03 | 保留原文 `hMxx/hTxx/hBxx/hIxx` 标签；marker 表见 04 号 §2.2 | 自动注释（SingleR/Azimuth）仅作交叉检查 | 原文 Fig. 1–2 | ★ |
| D9 | 恶性细胞 CNV | 03 | 按原文公式实现（100-gene moving average；baseline hE02–hE06；clip ±0.2）；`infercnv` 仅作敏感性 | 仅用 CNV score；仅用 CNV correlation | 原文 Methods「Identification of malignant cells」 | ★ |
| D10 | 轨迹算法 | 05 | 主路径 scVelo（velocyto）+ PAGA；URD 在 R 4.2 容器中复跑作敏感性 | 仅 URD；仅 PAGA | 原文 Methods；URD 版本兼容性 | ★ |
| D11 | 相关性网络基因集 | 06 | 用原文四条标准从 Smart-seq2 重新推导；若 Table S4A 可得则以其为准 | MSigDB 细胞类型 signature | 原文 Methods；Table S4A | ★ |
| D12 | LR 数据库与表达阈值 | 07 | FANTOM5 ∪ CellPhoneDB；>20% 细胞且 log2(TPM+1)>3；1,000 次置换 | 仅 CellPhoneDB；阈值 >10%/log2>2 | 原文 Methods | ★ |
| D13 | 小鼠–人类相似性 | 08 | 逻辑回归（elastic net α=0.99，10-fold CV，下采样至最小簇≥50）；人类 10x 为训练、小鼠为测试 | 互最近邻（Seurat `FindTransferAnchors`） | 原文 Methods「Similarity analysis」 | ★ |
| D14 | anti-CSF1R 丰度口径 | 09 | 仅 Study 2 用于簇频率定量（原文明确）；单侧 t 检验；单位 `% of CD45⁺` | 合并 Study1+2；双侧检验 | 原文 Methods；Fig. 4D 图注 | ★ |
| D15 | anti-CD40 时间点 | 10, 11 | day 2（早期激活）与 day 10（记忆/克隆）；tdLN 与肿瘤分开 | 仅 day 2；仅 day 10 | 原文 Fig. 5–7；S6C | ★ |
| D16 | TCR 与 STARTRAC | 11 | Cell Ranger 2.1.0 vdj；dominant αβ（UMI 最高）；UMI≥2；STARTRAC expansion/migration/transition | 仅克隆型描述 | 原文 Methods「TCR analysis」 | ★ |
| D17 | 生存分析分组 | 12 | TCGA COAD+READ；55th/45th 分位；Cox 校正 stage+MSI；HR+95% CI+n | 中位数二分；仅单因素 | 原文 Methods「Survival analysis」 | ★ |
| D18 | MSI 判定 | 12, 13 | 7 位点 PCR；≥2 单核苷酸不稳定 = MSI-H；公共数据用 TCGA 注释或 MSIsensor | 仅用 TCGA `MSI_status` | 原文 Methods；3 例 MSI-H | ★ |
| D19 | 湿实验范围 | 13 | S6 不执行湿实验；文档给 SOP 与验收标准；gap 表标"未执行" | 外包/合作实验室执行 | 本项目无样本与伦理批件 | ★ |
| D20 | 人–鼠 ortholog | 08 | 一标一 ortholog（biomaRt/OrthoDB，版本记录）；丢失基因列表落盘 | 直接用 symbol 交集 | 原文跨物种比较 | ★ |
| D21 | 多重比较口径 | 全部 | 组学默认 BH；原文只写 P<0.05 处双口径并列 | 仅 raw p | L0 §7 | ★ |
| D22 | GSVA 基因集 | 04, 09 | MSigDB C2 CP + Hallmark 主路径；GO BP 敏感性；原文未给版本 | 仅 Hallmark | 原文 Methods「GSVA」；图 3B 通路名 | ★ |

---

## 9. QC 门槛

- **数值层**：原文报告的关键数值与复现值允许 **±10%**（计数类 ±5%，p 值只比对数量级）；逐项见 14 号 §六。
- **判据层**：✅ 通过 → 推进；⚠️ 小幅偏差 → 排查后记录放行；❌ 显著偏差 → **第二套独立实现**复核，仍不符则入 gap 表。
- **结题层**：gap 表逐模块填「文献值 / 我方复现值 / 是否一致 / 差异原因 / 处理方式」。
- 全局 Go/No-Go（摘要）：
  - **G1** 细胞数：10x 43,817 ±5%；Smart-seq2 10,468 ±5%；anti-CD40 119,698 ±5%；anti-CSF1R 8,578 ±5%。
  - **G2** 簇数：10x 38 个免疫簇、Smart-seq2 36 个；13 个髓系、4 个 ILC、18 个 T、5 个 B；小鼠 15 个髓系簇。
  - **G3** hM12/hM13 marker 方向：C1QC↑/SPP1↓ 与 C1QC↓/SPP1↑ 正确；`TREM2`、`MERTK`、`CD80` 在 hM12 高；`MARCO`、`VEGFA` 在 hM13 高。
  - **G4** GSVA：SPP1⁺ 富集 angiogenesis/ECM；C1QC⁺ 富集 complement/antigen presentation（方向一致，FDR<0.05）。
  - **G5** 轨迹：FCN1⁺ → C1QC⁺/SPP1⁺ 方向；SPP1⁺ 与 NLRP3⁺ 相连。
  - **G6** 网络：TAM/cDC 节点度数最高；SPP1⁺–CAF/内皮、C1QC⁺–T 细胞边存在。
  - **G7** LR：CXCL10–CXCR3 在 C1QC⁺ 特异；SDC2–MMP2、SPP1–ITGAV、FN1–ITGA5 在 SPP1⁺ 特异（置换 p<0.05）。
  - **G8** 跨物种：mM15 ↔ hM13；mM11–14 ↔ hM12；mM06/07 ↔ hM04；mM04/05 ↔ hM03。
  - **G9** anti-CSF1R：mM12/mM14 显著减少；mM15 不减少；resistant 簇 `Vegfa`/`Cd274`/`Arg1` 升高。
  - **G10** 生存：C1QC^low SPP1^high 组 OS 最差（HR>1，p<0.05）；C1QC^high SPP1^low 组最好。
  - **G11** anti-CD40：mM07 cDC1-Ccl22 比例升高（day 2）；CD80/CD86/IL-12 升高；Bhlhe40⁺ Th1-like 与 CD8⁺ Tem 升高、Tex 降低。
  - **G12** 激活 DC signature 与 OS 正相关（HR<1，p<0.05；原文 HR=0.53，P=0.01）。
  - **G13** 阴性/不一致结果必须保留：未显著的簇、方向相反的亚组、无法复现的数值，均写入 gap 表。

---

## 10. 偏离声明（对 L0 的非红线偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| §8.1 头部四要素 | `00`/`14`/`15` 号文档不强制"模块对应原文"四要素（00 为总览、14 为质控、15 为背景） | L0 §1.3 自身将 00/NN-1/NN 定义为非模块文档 | 仅这三篇 |
| §6.4 下载校验 | 人类原始 FASTQ 不可得，只对 GEO processed matrix 与小鼠 processed matrix 做 md5 | 隐私限制；原文明确仅小鼠提供 FASTQ | `logs/00_checksums.md5` 覆盖范围在文件头注明 |
| §5.1 版本钉死 | 原文部分软件（URD 1.0.1、pagoda2 旧版）无法在 R 4.5 直接安装 | 版本断层；用容器/替代实现 | 05、02 号；D10、D5 |

> L0 §12 诚实性红线不可偏离。

---

## 11. 修订记录

| 版本 | 日期 | 修改内容 | 修改人 |
|---|---|---|---|
| v1.0 | 2026-10-18 | 首版：16 篇文档（00–15）结构、产物契约表、术语表、D 表索引（D1–D22）、环境版本表、全局 QC 门槛 G1–G13 | 复现项目组 |
