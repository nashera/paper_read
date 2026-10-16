# PROJECT_SPEC —— GLASS IDH 胶质瘤纵向演化 复现项目规范

> **上位规范**：`../REPRO_SPEC.md`（文献复现项目**统一规范（通用）v2.1**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（通用规范在本项目的实例化；不得违背 L0 红线）
> **本项目**：The GLASS Consortium. *Divergent genomic evolution in astrocytomas and oligodendrogliomas.* **Nature**（Open access），2026。PDF 标识 `s41586-026-11150-x`
> **数据**：Synapse **`syn17038081`**（Tables）+ **`syn26465623`**（Files）｜**代码**：**`github.com/verhaaklab/glass-idh`**（MIT）
> **项目目录**：`/mnt/e/paper_read/20261016_GLASS_IDH_glioma_divergent_evolution/`
> **启动日期**：2026-10-16　**本文件版本**：**v1.1**（对应 L0 v2.1）
>
> **效力**：L0 `../REPRO_SPEC.md` > **本文件 L1** > 13 号文档 §一 QC 门槛（L2）> 各模块文档（L3）。
> **D 表权威位置**：`13_复现质控清单、排错与批判性评述.md` **§七**（D1–D15）。

---

## 1. 项目定义

### 1.1 一句话目标

利用本文**公开的数据（Synapse）与公开的代码（GitHub）**，在 **L1 派生层**完整复现全部主图与统计量，并对不一致处如实报告；L2/L3（从 BAM/FASTQ 重跑）只做文档化。

### 1.2 三级复现层级（本项目一切取舍的依据）

| 层 | 输入 | 工具链 | 何时做 |
|---|---|---|---|
| **L1 派生层（主路径）** | Synapse 第 4 版衍生数据（变异/CN/临床） | `analysis/` + `figures/R/`（27 脚本）+ `tables/R/`（16 脚本）+ 我们的独立重算 | **默认全部在此层完成** |
| **L2 预处理层** | 受控访问的 BAM | `preprocessing_pipeline/`（251 文件；Snakemake 5.2.2 + GATK4 4.1.0.0 + SvABA + ASCAT…） | 仅当质疑官方变异/CN 检出时 |
| **L3 原始层** | FASTQ | 完整 GATK4 Best Practices | 一般不做 |

### 1.3 模块编号表

| 编号 | 文档 | 模块 | 层级 |
|---|---|---|---|
| 00 | 目录与复现总览 | 三级层级、索引、图表映射、全局门槛 | — |
| 01 | 数据与代码可得性、环境与仓库结构 | Synapse/SQL API、仓库结构、commit 锚定、环境契约、受控访问 | L1/L2/L3 |
| 02 | 队列构建、临床注释与样本QC | 581/239→231→183；QC1–QC6；barcode；183 vs 144；入组偏倚 | L1 |
| 03 | 变异检出与拷贝数调用 | GATK4/Mutect2/PoN/纵向 consensus；10-kb bin；GATK CNV + GISTIC 2.0.23 | L1/L2 |
| 04 | 驱动事件鉴定与互斥性分析 | dNdScv 三分集；GISTIC GeneGISTIC；DISCOVER；8 个扩增带 | L1 |
| 05 | 克隆性与纵向演化分析 | VAF→CCF；PyClone-VI；CCF>0.9；共享/空间共享比例 | L1/L2 |
| 06 | 突变负荷与超突变分类 | TMB 分母；超突变双条件；三态分类；MSH6/MMR；TMZ | L1 |
| 07 | 突变签名分析与治疗归因 | NMF；COSMIC v3.5 + Signal；Extractor/Assignment；Palimpsest | L1 |
| 08 | 拷贝数签名与全基因组加倍 | ASCAT + Sequenza 双过滤；48-feature CN 矩阵；SigProfilerAssignment | L1/L2 |
| 09 | 非整倍体与FGA量化 | 臂级 80% 规则；段级 ±0.2；FGA 公式与平台分母 | L1 |
| 10 | 染色体不稳定性 | SvABA；ShatterSeek；kataegis；AmpliconSuite ecDNA | L1/L2 |
| 11 | 生存分析与预后模型 | 四终点；KM/log-rank；Cox+Schoenfeld；Harrell C+bootstrap | L1 |
| 12 | 外部验证与批次效应检验 | TCGA-LGG/GENIE v18.5/CBTRUS；分子替代分型；分层 logistic | L1 |
| 13 | 复现质控清单、排错与批判性评述 | QC 表、排错、原文 25 条疏漏、批判性评述、**D 表**、gap 表 | — |
| 14 | 背景知识精讲与方法学原理速查 | IDH 生物学、CIN 现象学、签名体系、克隆演化、工具速查 | — |

---

## 2. 目录与命名

### 2.1 目录

```
raw/{synapse,glass_raw,reference,cosmic,genie}/   # 原始与参考数据（只读）
processed/     # 中间对象（tsv/rds/parquet/npz）
results/       # 表格结果（csv/tsv/rds）
figures/       # 出图（pdf + png，dpi=300）
code/          # 我们的编排与改造脚本
code/glass-idh/  # ★ 官方仓库 clone（只读、不改、记录 commit hash）
logs/          # DECISIONS.md / sessionInfo / 校验和 / Synapse 查询快照
```

### 2.2 命名特例

| 情形 | 规则 |
|---|---|
| 原文 PDF | 保留原题全称（含空格），冒号改连字符 |
| 图表 | 与原文图号对齐：`Fig<N><L>_<语义>.pdf`、`ED_Fig<N><L>_<语义>.pdf`（如 `ED_Fig3_forest_dnds_ci.pdf`、`Fig3e_prs_cdkn2a_focalamp.pdf`） |
| **核心共享对象** | **无编号前缀 + `glass_` 语义前缀**（如 `processed/glass_pgeno.tsv`、`glass_cnv_calls.tsv`、`glass_gold_pairs.tsv`）—— 被 ≥2 模块消费 |
| 模块私有产物 | `<NN>_` 前缀（`results/09_fga.csv`、`results/11_km_logrank.csv`） |
| 官方脚本 | 引用时用 `figures/R/<脚本名>`；**不确定的写"以仓库当前版本为准"**，不得凭记忆编造 |
| 可复用脚本 | 放仓库 `../tools/`（不复制进项目） |

### 2.3 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `data/processed/`、`data/raw/` | `processed/`、`raw/` | L0 §3.1 |
| `processed/glass_geno.tsv`（若出现） | `processed/glass_pgeno.tsv`（患者级基因型矩阵） | L1 本文件 |
| `results/09_aneuploidy_*.csv` 与 `09_arm_level_calls.csv` 并存 | 统一为 `09_arm_level_calls.csv`（臂级）、`09_segment_level_counts.csv`（段级）、`09_fga.csv`（全基因组） | L1 本文件 |
| 把 `legacy_codes/` 里的脚本名当现行脚本 | **一律作废**；现行脚本只在 `analysis/`、`figures/`、`tables/` | 官方声明 |
| `processed/cohort_183_pairs.csv`（18 处引用） | ✅ **`processed/glass_gold_pairs.tsv`**（L1 权威名，供分析代码使用）；CSV 版降级为**同生产者导出的别名** | L1 本文件（消解 02 号与 03–11 号之间的契约冲突） |
| `processed/02_paired_manifest.tsv` | 02 号**模块私有中间件**（非跨模块契约），不得被其他模块消费 | L1 本文件 |

---

## 3. 产物契约表

> 纪律（L0 §4）：**每个产物唯一生产者**；改格式/名 = 改本表 + 通知消费方 + 记 `logs/DECISIONS.md`。

### 3.1 核心共享对象（`glass_` 前缀，跨模块）

| 产物 | 生产者 | 格式 / 结构（下游必须知道的） | 消费方 |
|---|---|---|---|
| `processed/glass_gold_pairs.tsv` | 02（**唯一生产者**） | **金标准配对表（L1 权威名）**。必需列：`case_barcode, tumor_pair_barcode(29 字符), tumor_barcode_a, tumor_barcode_b`（后三者**不在 Supplementary Table 1 中，须由 02 依据 barcode 规则合成**）、`surgery_index_initial, surgery_index_recurrent, interval_days, type(astro/oligo), platform(WGS/WXS), qc_flag` | 03–12 |
| `processed/cohort_183_pairs.csv`（**别名**） | 02（同一生产者） | 与上表**同内容的 CSV 导出**，仅供人工查看与外部交接 | 人工/交接 |
| `processed/glass_surgeries.tsv` | 02 | 手术级临床表（含治疗、分级、日期） | 02/06/11/12 |
| `processed/glass_pgeno.tsv` | 03 | 患者级基因型/驱动矩阵（基因 × 患者，值 = 事件类型/时间点） | 04/06/08/11 |
| `processed/glass_cnv_calls.tsv` | 03 | 基因级拷贝数调用（含 `log2_ratio`、`call(amp/hd/neutral)`、权重均值/±2SD 阈值记录） | 04/08/09 |
| `processed/glass_gatk_seg.tsv` | 03 | GATK 分割段（`chr, start, end, log2_ratio, sample`）—— **b37/hg19，染色体无 `chr` 前缀或统一规则必须记录** | 08/09/10 |
| `processed/glass_seg_stats.tsv` | 03/04 | 每样本的段数、中性段比例等派生于分割的统计 | 09 |
| `processed/glass_hm_status.tsv` | 06 | 超突变**三态**状态（`HM`/`non-HM`/`unassigned`）+ SBS11 比例 | 07/11/12 |
| `processed/05_mutation_pair_timing.tsv` | 05 | 突变级时序（initial-only / recurrent-only / shared）+ VAF + CCF | 05/07/11 |

### 3.2 模块产物索引

| 模块 | 关键产物 |
|---|---|
| 01 | `logs/01_synapse_query_snapshot.tsv`（表名/行数/查询时间）、`logs/01_repo_commit.txt`（`git rev-parse HEAD`）、`logs/01_sessionInfo.txt` |
| 02 | `results/02_baseline_table.csv`、`02_pair_structure.csv`、`02_inclusion_bias.csv`、`02_batch_effect.csv` |
| 03 | `processed/glass_*`（见 3.1）、`results/03_qc_coverage.csv`、`03_platform_strata.csv` |
| 04 | `processed/driver_changes.tsv`、`results/04_dndscv_by_stratum.csv`、`04_gistic_focal_events.csv`、`04_discover_mutex.csv` |
| 05 | `processed/05_pyclone_input.tsv`、`05_crosscounts.tsv`、`results/05_acquired_driver_clonality.csv`、`05_ccf_disagreement.csv` |
| 06 | `results/06_hypermutation_status.tsv`（+ `__strict` 敏感性）、`06_ed8a_tmb_grade.csv`、`06_fig3g_paired_tests.csv`、`06_ed8d_interval_correlations.csv` |
| 07 | `results/07_crosscheck_cos_sim.csv`、`07_excluded_4_cases.csv`、`07_paired_wilcoxon_initial_vs_recurrence.csv`、`07_tmz_dose_response.csv` |
| 08 | `results/08_wgd_sequenza.csv`、`08_wgd_final_by_timepoint.csv`、`08_wgd_mcnemar.csv`、`08_cnsig_relative_activity.csv`、`08_cnsig_by_wgd.csv` |
| 09 | `results/09_arm_level_calls.csv`、`09_fga.csv`、`09_fga_paired_tests.csv`、`09_aneuploidy_score.csv`、`09_sensitivity.csv` |
| 10 | `results/10_chromothripsis_calls.csv`、`10_chromothripsis_approx.csv`、`10_kataegis_approx.csv`、`10_ecdna_heuristic.csv`、`10_cin_matrix.csv` |
| 11 | `processed/11_patient_level.tsv`、`11_paired_endpoints.tsv`、`results/11_km_logrank.csv`、`11_cox_uni.csv`、`11_cox_multi.csv`、`11_cindex_bootstrap.csv`、`11_baseline.csv` |
| 12 | `processed/12_genie_tidy.rds`、`12_tcga_tidy.rds`、`12_unified_cohort.tsv`、`results/12_batch_logistic.csv`、`12_cbtrus_table.csv` |
| 13 | `logs/DECISIONS.md`、gap 表 |

### 3.3 官方脚本对照表（**2026-10-16 经 GitHub API 实测枚举，可直接引用**）

> 下列名字来自 `https://api.github.com/repos/verhaaklab/glass-idh/contents/<path>` 的实际返回。
> **脚本名确凿**；但"脚本 ↔ 具体图号"的对应为**推断**（脚本头部注释未逐一核对），引用时以仓库中的脚本标题注释为准。

**`figures/R/`（27 个 R 脚本）**

| 推断对应的原文图 | 脚本名 |
|---|---|
| Fig. 1（队列/时间线/治疗） | `GLASS-I_SwimmersPlot.R`、`GLASS-I_Cohort_Description.R`、`GLASS-I_TMZ_Cycles.R`、`GLASS-I_IDH_Mutations.R` |
| ED Fig. 1（QC 流程） | `GLASS-I_Sankey_QC.R` |
| Fig. 2（驱动景观） | `GLASS-I_dNdScv.R`、`GLASS-I_sSNV_Mutual_Exclusivity.R`、`GLASS-I_CN_Heatmap.R` |
| Fig. 3（生存） | `GLASS-I_Univariable_KM_Hypermutation.R`、`GLASS-I_Univariable_KM_sSNV.R`、`GLASS-I_Univariable_KM_CNA.R`、`GLASS-I_Univariable_KM_Drivers_GENIE.R`、`GLASS-I_Mutation_Burden_NHM.R`、`GLASS-I_Univariable_KM_ID8.R` |
| Fig. 4（签名） | `GLASS-I_Mutational_Signature_Boxplots.R`、`GLASS-I_TMZ_Cycles.R` |
| Fig. 5（CIN / 拷贝数） | `GLASS-I_FGA.R`、`GLASS-I_FGA_NHM.R`、`GLASS-I_FGA_NHM_Grade.R`、`GLASS-I_CN_Segments_{Astro,Oligo}_{Primary,Recurrence,Delta}.R`（6 个）、`GLASS-I_CIN_Markers.R`、`script_landscape_chromosomal_instability.R` |
| ED Fig. 2（与 TCGA 比较） | `GLASS-I_TCGA_OS_Comparison.R` |

（未列入上表的脚本：`GLASS-I_CN_Segments_*.R` 已列；此外 `figures/R/` 内含 `.gitkeep`。）

**`analysis/R/`（12 个）**：`GLASS-I_Amplicon_Architect.R`、`GLASS-I_Chromothripsis.R`、`GLASS-I_ShatterSeek.R`、`WGD_sequenza_script.R`、`ascat_analysis.R`、`ascat_merge.R`、`build_idh_mut_gold_wgs_tel_length_norm.R`、`cnsig_expl.R`、`create_idh_mut_gold_purity_coverage_from_files.R`、`prepare_gistic_input.R`、`script_fga_calculation.R`、`sigp_matgen_cnv_prepare.R`

**`tables/R/`（16 个）**：`GLASS-I_Baseline.R`、`GLASS-I_Drivers_Prior_Treatment.R`、`GLASS-I_Drivers_Table.R`、`GLASS-I_IDH_Mutations_Table.R`、`GLASS-I_McNemar_CNA.R`、`GLASS-I_Multivariable_{OS,PRS,TTR}.R`、`GLASS-I_Prior_Treatment_Table.R`、`GLASS-I_Summarizing_Table_Molecular_Cohort.R`、`GLASS-I_Univariable_{OS,PRS,TTR}_{Clinical,Molecular}.R`

**`analysis/SQL/`（4 个）**：`GLASS-I_CNA_Data.sql`、`GLASS-I_dNdScv_HM_input.sql`、`GLASS-I_dNdScv_NHM_input.sql`、`gistic_prepare.sql`

**`figures/python/`（1 个）**：`GLASS-I_World_Map_Cohort_Distribution.py`

**其他确凿路径**：`analysis/bash/mutational_signature_pipeline.sh`、`analysis/bash/kataegis_pipeline.sh`、`analysis/python/ampsuite_manifest.py`、`analysis/slurm/GLASS-I_SvABA.sbatch`、`analysis/snakemake/{ampsuite,ascat,sigprofilertoolkit}.smk`

> ⚠️ **规范要求**：模块文档引用官方脚本时**必须使用本表的名字**；若本表未覆盖，写"以仓库当前版本为准（`ls figures/R/` 枚举）"，**禁止凭记忆编造脚本名**。


---

## 4. 术语与判定口径（唯一化）

| 术语 | 定义（**唯一口径**） |
|---|---|
| 分子队列 | **183 例 = 94 astrocytoma + 89 oligodendroglioma**（QC 后，两时间点） |
| 临床注释集 | **119 例 astro（262 手术）/ 112 例 oligo（292 手术）** —— **与"分子队列"不是同一数字**，禁止混用 |
| 配对口径 | **`initial–recurrence`（n=183）** vs **`primary–recurrent`（n=144）** —— 每张图用哪套必须核图注（D5） |
| 平台 | `WGS`（39 astro / 56 oligo）与 `WXS`（55 astro / 33 oligo）**分平台报告**，FGA/TMB 不可混算 |
| HM | `HM`（>10 mut/Mb 且 SBS11 >25% 复发特异突变）/ `non-HM`（<10 mut/Mb）/ **`unassigned`（>10 但 SBS11 <25%）** |
| clonal | PyClone-VI 的 **CCF > 0.9** |
| FGA | $FGA = 1 - \frac{\text{neutral bp}}{\text{callable bp}}$；**WGS 分母 3×10⁹ bp、WXS 5×10⁷ bp** |
| 臂级非整倍体 | 合并段覆盖 **>80%** 臂长才给 −1/0/+1；**±0.1** 判增减、**±0.05** 判变化 |
| 段级非整倍体 | **\|log2 ratio\| > 0.2** 记为 altered |
| WGD | **ASCAT 多样本判定** 且 **Sequenza major CN fraction ≥ 0.5**（双重过滤） |
| 参考基因组 | **b37/hg19（human_g1k_v37_decoy）**；SvABA 用 **hs37d5**；**不切 hg38**（D4） |
| OS / TTR / OS-from-birth / PRS | 四个终点定义见 11 号文档；**PRS 起点 = QC 通过的末次复发手术** |
| GENIE 分型 | **IDH 突变且无 TP53/ATRX 改变 = oligodendroglioma**；含改变者 = astrocytoma |
| **配对口径 `144`** | **`primary–recurrent` = 配对表中 initial 手术号 == 1**（=122+20+2）—— 原文只说"indicated accordingly"，**未给定义**；由 Supplementary Table 1 实测确认（D5） |
| **case barcode 前缀** | **不总是 `GLSS`**：实测 183 例中 `GLSS` 88、`OLIG` 73、`CARE` 10、`TCGA` 12 → **`grepl("^GLSS-")` 会静默丢 95 例**；`tumor_pair_barcode` 为 **29 字符**（非 35） |
| **Data Release 混淆** | **GLASS Data Release 与肿瘤类型/平台强相关**（实测 2022：65 astro/29 oligo 且 30 WGS/64 WXS；2026：29 astro/60 oligo 且 65 WGS/24 WXS）→ 批次检验**必须**按类型+超突变分层，否则 release/platform 效应会混入"中心效应"（D12） |
| **数据可达性** | **本机对 Synapse 被地理封锁**（GEO_RESTRICTION / 403），`genie.cbioportal.org` 451；**公开替代通道 = Nature 补充材料 `41586_2026_11150_MOESM1_ESM.xlsx`（114,368 B，sha256 `13b4e591…`，12 sheets = Supp. Tables 1–12，无需账号）**（D22） |

---

## 5. 环境与版本钉死

| 组件 | 版本 | 说明 |
|---|---|---|
| R | **4.3.2** | 原文 |
| Python | **3.10.13** | 原文 |
| 流程 | Snakemake **5.2.2**；GATK4 **4.1.0.0**；BCFtools **1.9** | L2 用 |
| 关键工具 | GISTIC **2.0.23**；dNdScv **0.1.0**；SvABA **1.2.0**；ShatterSeek **1.1**；SigProfilerExtractor **1.2.1**；SigProfilerClusters **1.2.2**；Palimpsest **2.0.0**；DISCOVER **3.1.7**；PyClone-VI；ASCAT；Sequenza；AmpliconSuite（AmpliconArchitect + AmpliconClassifier） | 版本偏离须记 `DECISIONS.md` |
| 参考库 | **COSMIC v3.5**（SBS1–96 + 全部 ID + CN 签名）+ **Signal**（SBS97+） | COSMIC 需注册/许可 |
| 数据版本 | **Synapse 第 4 次发布**（记录查询日期 + 表名 + 行数） | D3 |

| 代码版本 | `code/glass-idh/` 的 **commit hash** —— 实测 **HEAD = `6bdd7292`（2026-09-18）**；**经核实该仓库不是 `fsvarn/GLASSx` 的 fork**（694 commits、起始 2018-05-01、`.git/config` 无 upstream；README 把 GLASSx 列为既往论文代码而已） | D2 |

> ⚠️ **溯源声明（PI 追加，2026-10-16）**：`syn17038081`、`syn26465623`、`github.com/verhaaklab/glass-idh` 三者的**出处均为原文的 Data availability / Code availability 段**（前者原文："Deidentified somatic variant and clinical datasets are publicly available through Synapse with Project SynID syn17038081"；后者给出 GitHub 地址）。
> - **GitHub 仓库：已由本机通过 GitHub API 实测核实**（`figures/R` 27 个脚本、`analysis/R` 12 个、`tables/R` 16 个的真实文件名见 §3.3，均为实际返回）。
> - **Synapse：PI 未能独立核实** —— 本机对 `repo-prod.prod.synapse.org` 的 DNS 解析失败（ENOTFOUND），网页端返回 **403**（需登录/反爬拦截）。
>
> **因此，正式开工前使用者必须在浏览器登录 Synapse 后人工确认三件事**：① 项目 `syn17038081` 存在且当前账号可读；② 其下 `tables/` 的**实际表名与行数**；③ 是否需签署数据使用协议（DUA）。
> 在上述确认完成前：**不得**在文档中写入任何具体表名，**不得**把"数据已就绪"写进任何报告（L0 §12 诚实性红线）。核实结果记入 `logs/DECISIONS.md`。


---

## 6. 统计与随机性

| 项 | 取值 |
|---|---|
| 随机种子 | R `set.seed(20260101)`；Python `np.random.seed(20260101)`（NMF replicates、bootstrap、模拟基因组、PyClone） |
| 非配对分类 | Fisher 精确 / χ² |
| 非配对连续 | Wilcoxon rank-sum / Kruskal–Wallis |
| **配对**分类 | **McNemar** |
| **配对**连续 | **Wilcoxon signed-rank** |
| 生存 | KM + log-rank；Cox + **Schoenfeld** 验 PH；似然比 P≤0.10 纳变量 + 逐步后向；**Harrell C + 1000 bootstrap** |
| 多重校正 | 关键结论 **BH(FDR)**；与原文 Supp. Table 10 口径对齐（D15） |
| 显著性 | 双侧 **P < 0.05** |
| 报告完整性 | 每个比较必须给出 **n、检验名、效应量/HR 与 95% CI**；**方向与原文不一致必须并列报告** |

---

## 7. 决策登记表（D 表）索引

> **权威位置 = `13_复现质控清单、排错与批判性评述.md` §七**（**D1–D24**）。本表仅列 ★（改变结论口径）项；**D16–D24**（生存亚组 n 口径、签名笔误、治疗分母、MSH6 计数、Table 1 随访口径、治疗暴露口径、**Synapse 地理封锁应对**、**参考库三套版本**、**公开代码链路缺口**）见该处。

| # | 决策项 | 一句话 | ★ |
|---|---|---|---|
| D1 | 复现层级 | 只做 L1；L2/L3 文档化 | ★ |
| D2 | 代码锚定 | 记录 commit hash；核 fork 与论文版本一致性 | ★ |
| D3 | 数据版本锁定 | 锁第 4 次发布；表名以页面为准 | ★ |
| D4 | 参考基因组 | 保持 b37/hs37d5，不切 hg38 | ★ |
| D5 | 配对口径 | 按图注选 183/144，并做"仅 Tx–R1"敏感性 | ★ |
| D6 | FGA 平台口径 | WGS/WXS 严格分层 | ★ |
| D8 | TMB 分母与过滤 | 分平台；覆盖 ≥15× | ★ |
| D9 | 超突变三态处置 | unassigned 单列，不并入 non-HM | ★ |
| D10 | 签名参考库 | COSMIC v3.5 + Signal；仅 WGS | ★ |
| D12 | GENIE 分子替代分型 | 加错分模拟 | ★ |
| D13 | 生存模型稳健化 | 原文逐步法 + bootstrap/惩罚 Cox | ★ |
| D14 | PRS 条件性选择 | 从 R1 起算的敏感性分析 | ★ |
| D15 | 多重比较 | 关键结论补 BH 并报告 | ★ |
| D25 | 配对表契约收敛 | 分析代码一律用 `processed/glass_gold_pairs.tsv`；CSV 版为别名 | ★ |

---

## 8. QC 门槛

- **全局门槛**：`00_目录与复现总览.md` §五（9 条，含 183/94/89、122+61、7.8/21.4 年、38/94、16/70 与 20/74、Fig. 4c 六个 P 值、WGD McNemar 0.55/0.006 等）。
- **逐项门槛**：`13_...md` §一（**45 项** QC 表，含层级归属与判据）。
- **L1 最强验收**：从 Synapse 公开数据**独立重算**全部统计量，p 值比量级、比例 ±10%、**方向必须一致**。
- **结题**：`13_...md` §六 gap 表。

---

## 9. 偏离声明（对 L0 非红线规则的偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| L0 §3.2 禁空格/中文 | 原文 PDF 与 `NN_*.md` 保留空格与中文 | L0 明文豁免 | 无 |
| L0 §10.1 D 表位置 | D 表置于 13 号文档 | L0 允许"二选一"；13 号已被大量引用 | 需同步本文件 §7 |
| L0 §3.2 前缀规则 | 核心共享对象用 **`glass_` 语义前缀**而非无前缀 | 与官方仓库的命名风格一致、便于与 `NN_` 模块产物区分 | 已在 §2.2/§2.3 规范化 |

---

## 10. 修订记录

| 版本 | 日期 | 修订内容 | 说明 |
|---|---|---|---|
| v1.1 | 2026-10-16 | 按审核结果修订：① 纠正"`glass-idh` 是 `fsvarn/GLASSx` 的 fork"（实为独立仓库；记录 HEAD `6bdd7292`）；② §5 溯源声明移到表格之后（修复表格被劈开的渲染 bug）；③ §3.3 补 `analysis/SQL`（4）与 `figures/python`（1）；④ §4 术语表新增 `144` 定义、**barcode 前缀不总是 GLSS**、**Data Release 与类型/平台强相关**、**Synapse 地理封锁与公开替代通道**；⑤ §7 D 表索引扩到 D1–D24（新增地理封锁应对、参考库三套版本、公开代码链路缺口等） | 本版本 |
| v1.0 | 2026-10-16 | 首版：三级复现层级（L1 主路径）、模块编号表、**核心共享对象契约（`glass_` 前缀）**、术语与判定口径唯一化（含"183 vs 119/112""183 vs 144""WGS/WXS 分层"）、环境与数据/代码版本锚定、统计约定、D 表索引（D1–D15）、L1 最强验收判据、偏离声明 | 本版本 |
