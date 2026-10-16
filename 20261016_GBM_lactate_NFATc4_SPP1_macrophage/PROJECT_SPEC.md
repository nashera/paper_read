# PROJECT_SPEC —— GBM 乳酸代谢 / NFATc4–SPP1 轴 复现项目规范

> **上位规范**：`../REPRO_SPEC.md`（文献复现项目**统一规范（通用）v2.1**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（通用规范在本项目的实例化；不得违背 L0 的红线，可填写 L0 未规定的项目特有内容）
> **本项目**：Wang R, Wu W, Wu X, Yuan Y, Zhang Z, Chen J, Xie Y, Zhang L, Yu X, Zhang H, Xie Z. *Integrative single-cell and spatial transcriptomic analysis reveals a lactate-driven crosstalk between NFATc4⁺ tumor cells and SPP1⁺ macrophages in glioblastoma.* **J Transl Med (2026) 24:1253.** DOI: [10.1186/s12967-026-08555-7](https://doi.org/10.1186/s12967-026-08555-7)
> **项目目录**：`/mnt/e/paper_read/20261016_GBM_lactate_NFATc4_SPP1_macrophage/`
> **启动日期**：2026-10-16　**本文件版本**：**v1.1**（对应 L0 v2.1）
>
> **效力**：本项目内，L0 `../REPRO_SPEC.md` > **本文件 L1** > 13 号文档 §一 QC 门槛（L2）> 各模块文档（L3）。
> 任何模块文档与本文件冲突时，以本文件为准；本文件与 L0 红线冲突时，以 L0 为准。

---

## 1. 项目定义

### 1.1 一句话目标

在公开数据上**完整、可核验地重建**原文证据链：单细胞乳酸代谢异质性 → pySCENIC 锁定 NFATc4 → 空间共定位 SPP1⁺ 巨噬细胞 → 六队列生存 → 101 组合机器学习预后模型 → 免疫/基因组/药敏延伸 → 体外功能验证；并对原文未写清处给出显式复现决策（D 表）。

### 1.2 科学问题（原文答案链，复现时逐条验证）

1. GSE138794 单细胞图谱：**12,618 细胞 / 24 cluster / 10 大类细胞**，肿瘤细胞与巨噬细胞互作最强（Fig. 1–2）。
2. **332 LRGs + ssGSEA** 划分高/低乳酸肿瘤细胞；pySCENIC 得 **6 个 regulon 模块**，高乳酸组升高 TF：**SPI1、NFATc4、HOXC9、POU3F3、FOXO1**（Fig. 3, S1）。
3. Monocle 拟时序：**2 个分支点 / 5 个 state**，NFATc4 沿第二分支升高；**AHCYL1 先升后降、CADPS 持续上调**（Fig. S2）。
4. 空间转录组 GSE194329（**4 切片 / 11,970 spot**）：NFATc4⁺ 肿瘤细胞与 SPP1⁺ 巨噬细胞显著共定位，高乳酸区距离更近（Fig. 4, S3）。
5. 六队列中高 NFATc4 普遍预后不良；Mayo PDX 中两个 signature 强正相关（Fig. 3D, S4）。
6. 体外：Nfatc4 敲低抑制 GL261 增殖/迁移、诱导凋亡与周期阻滞；Spp1 敲低在巨噬细胞中同样抑制共培养肿瘤细胞；乳酸上调 NFATc4（Fig. 5, 6）。
7. 机器学习：**StepCox[forward] + RSF** 最优（TCGA C-index **0.87**、GSE4412 **0.64**、meta **0.62**；5 年 AUC **0.976/0.905**；pooled HR **1.82, 95% CI 1.64–2.02**）（Fig. 7–9, S5）。

### 1.3 交付物清单

| 类别 | 文件 | 状态 |
|---|---|---|
| 原文 | `Integrative single-cell and spatial transcriptomic analysis ... (J Transl Med 2026, 24-1253).pdf` | ✅ 已归档 |
| 项目规范 | **`PROJECT_SPEC.md`（本文件，L1）** | ✅ v1.0 |
| 通用规范 | `../REPRO_SPEC.md`（L0） | ✅ v2.1 |
| 模块文档 | `00`–`14` 共 15 篇（约 1.2 MB / 15,300 行） | ✅ 已产出并交叉审核 |
| 代码 | `code/*.{R,py,sh}`（编号前缀与文档一致） | ⏳ 待实现 |
| 产物 | `raw/ processed/ results/ figures/ logs/` | ⏳ 待生成 |

### 1.4 模块编号表

| 编号 | 文档 | 模块 |
|---|---|---|
| 00 | 目录与复现总览 | 索引、流程图、数据资产、全局门槛 |
| 01 | 环境搭建与数据获取 | 双环境、版本表、6 类数据源、**332 LRG 基因集重建** |
| 02 | 单细胞QC整合与十大类细胞注释 | QC、DoubletFinder、Harmony、10 大类、inferCNV、CytoTRACE |
| 03 | 乳酸代谢评分与高低乳酸亚群划分 | 332 LRGs + ssGSEA/AUCell、分组阈值、与 glycolysis/hypoxia 解耦 |
| 04 | pySCENIC 与 NFATc4 发现 | GRNBoost2 → RcisTarget → AUCell、6 模块、差异 TF |
| 05 | 拟时序轨迹与分支依赖基因 | Monocle、2 分支 5 state、BEAM、GO |
| 06 | 空间反卷积与共定位 | SCTransform → RCTD、ULM、squidpy 邻域富集、距离统计 |
| 07 | 细胞通讯与 NFATc4–SPP1 轴 | CellChat（scRNA + 空间）、SPP1 通路、循环论证的正交敏感性分析 |
| 08 | bulk 队列整合与 NFATc4 生存分析 | 七队列构建、KM/Cox、截断点偏差、Mayo PDX |
| 09 | 机器学习 101 组合预后模型 | 共识聚类 → DEG → unicox → Mime → StepCox[forward]+RSF |
| 10 | 免疫浸润、突变与 CNV | IOBR、maftools、GISTIC 2.0、GSVA |
| 11 | 药物敏感性与 CMAP 重定位 | CMAP_XSum、oncoPredict、CTRP/PRISM、scMetabolism |
| 12 | 体外功能验证实验方案 | siRNA、RT-qPCR、WB、EdU/凋亡/周期/Transwell、共培养 |
| 13 | 复现质控清单、排错与批判性评述 | QC 表、排错、**D 表（本项目权威位置）**、批判性评述 |
| 14 | 背景知识精讲与方法学原理速查 | 生物学/算法原理、术语、46 篇文献导读 |

---

## 2. 目录与命名（本项目的具体化）

### 2.1 数据源细分

```
raw/
├── geo/{GSE138794,GSE194329,GSE108474,GSE4412}/
├── tcga/            # star_counts / tpm / clinical
├── cgga/            # mRNAseq_693 / mRNAseq_325 / mRNA-array_301
├── cbioportal/      # Mayo Clinic GBM PDX
├── cistarget/       # motif 注释 feather + tbl（物种/基因组版本必须匹配）
├── spatial/GSE194329/<GSM>/   # Space Ranger outs/（含 tissue_positions 与 scalefactors）
└── drug_response/   # GDSC / CTRP / PRISM / LINCS
```

### 2.2 项目命名特例

| 情形 | 规则 |
|---|---|
| 原文 PDF | 保留原题全称，**冒号改连字符**（`24:1253` → `24-1253`），允许空格 |
| 图表 | 与原文图号严格对齐：`Fig<N><L>_<语义>.pdf/png`（如 `Fig3D_NFATc4_KM_six_cohorts.pdf`、`FigS4B_pdx_signature_scores.csv`） |
| 队列占位 | `<cohort> ∈ {TCGA, CGGA693, CGGA325, CGGA301, GSE108474, GSE4412, meta}`（**CGGA325 替代原文 CGGA625**，见 D5） |
| 切片占位 | `<GSM> ∈ {GSM5833533, GSM5833534, GSM5833535, GSM5833536}` |
| 模块私有产物 | `<NN>_` 前缀（`02_celltype_proportion.csv`、`06_rctd_weights_<GSM>.csv`） |

### 2.3 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `data/raw/`、`data/processed/` | `raw/`、`processed/` | D2 |
| `processed/seurat_gbm_annotated.rds` | `processed/seurat_gbm_qc_seurat5.rds` | D1 |
| `processed/bulk_CGGA.rds` | `processed/bulk_CGGA693.rds` / `_CGGA325.rds` / `_CGGA301.rds` | D5 |
| `results/model_coef.rds`、`results/model_final_coefficients.csv` | `results/ml/stepcox_forward_coefficients.csv` | D4 |
| `results/risk_score_tcga.csv` | `results/ml/risk_score_all_cohorts.csv` | D4 |
| `figures/Fig3D_KM_six_cohorts.pdf` | `figures/Fig3D_NFATc4_KM_six_cohorts.pdf` | L1 本文件 |
| `results/05_BEAM_branch2_qval<0.01.csv` | `results/05_BEAM_branch2_qval_lt_0.01.csv` | L1（禁 `<`） |
| 正文/表格中的 `..._zscore_{s}.csv` | `..._zscore_<GSM>.csv`（**仅文档表述**；代码 f-string 保持 `{s}`） | L1 |
| 原文 "CGGA625" | `CGGA325 (CGGA625 replacement)` | D5 |

---

## 3. 产物契约表

> 纪律：**每个产物只有一个生产者**（L0 §4）。改格式/名 = 改本表 + 通知全部消费方 + 记 `logs/DECISIONS.md`。

### 3.1 核心共享对象（无前缀）

| 产物 | 生产者 | 格式 / 结构（下游必须知道的一切） | 消费方 | 脚本 |
|---|---|---|---|---|
| `processed/seurat_gbm_qc_seurat5.rds` | 02 | Seurat 5 对象；`celltype`（10 level，见 §4）、`seurat_clusters`（24）、`sample`（5 GSM）、`cnv_score`/`cnv_call`、`S.Score`/`G2M.Score`、`cytotrace`；（v3 Assay 版为 `seurat_gbm_qc.rds`） | 03/04/05/06/07 | `code/02_*.R` |
| `processed/02_cell_metadata.tsv` | 02 | 同上元数据的轻量表（barcode 为行名） | 03–07（快速读取替代） | `code/02_*.R` |
| `processed/lrgs_332.txt` | 01 | `character(332)`，HGNC symbol，**txt 非 rds**；来源分列见 `processed/lrgs_source_lists.rds` | 03/04/08/09/14 | `code/01_build_LRGs.R` |
| `processed/neoplastic_lactate_scored.rds` | 03 | **已 subset 到 neoplastic**；`lactate_score`（主，ssGSEA）、`lactate_group`（factor，levels `c("Low","High")`，**中位数二分**）、`lactate_aucell`、`lactate_ucell` | 04/05/07 | `code/03_*.R` |
| `processed/pyscenic_input.loom` | 04 | TF ∪ HVG 表达矩阵（log 化）；细胞与 03 的 neoplastic 子集一致 | 04 | `code/04_01_*.R/py` |
| `results/04_aucell_mtx.csv` | 04 | cell × regulon 的 AUCell 活性；**列名带 `(+)` 后缀**；需自行投影到 Seurat 对象 | 05/07 | `code/04_03_*.R` |
| `processed/06_visium_sct_list.rds` / `06_visium_annotated_list.rds` | 06 | 4 切片 Visium 对象（list，按 GSM 命名）；含 RCTD 权重、区域标签（high-lactate/low-lactate/non-tumor） | 07 | `code/06_*.R` |
| `processed/bulk_matrix.rds` / `processed/bulk_<cohort>.rds` | 01 | `list(exp = 基因×样本 matrix, surv = data.frame(OS.time, OS), anno)`；**OS.time 统一为天**；exp 未 log2 | 04/08/09/11 | `code/01_*.R` |
| `processed/meta_cohort.rds` | 01 | 合并队列，**构成 = M1（不含 TCGA，见 D3）**；含 `combat` 版本矩阵 | 09 | `code/01_meta_cohort.R` |
| `processed/pdx_mayo.rds` | 01 | Mayo Clinic GBM PDX 表达 + 临床 | 08 | `code/01_pdx_cbioportal.R` |
| `processed/sig_spp1_mac_tam.rds` / `sig_nfatc4_tumor.rds` | 07 | `FindMarkers` 得到的 signature 基因表；**已剔除 SPP1 / NFATC4 自身**（避免自证） | 08（PDX ssGSEA） | `code/07_*.R` |
| `results/ml/stepcox_forward_coefficients.csv` | 09 | 列 `gene, coef, HR` | 11/14 | `code/09_*.R` |
| `results/ml/risk_score_all_cohorts.csv` | 09 | 列 `sample, cohort, riskScore, group` | 11/14 | `code/09_*.R` |

### 3.2 模块产物索引

| 模块 | 主要产物（节选；完整清单见各文档 ⑧ 图表对应表） |
|---|---|
| 01 | `raw/**` 全量 + `logs/geo_manifest.tsv` + `logs/00_checksums.md5` + `logs/01_lrg_membership.tsv`（NFATC4 循环论证检查，D11） |
| 02 | `02_qc_pre/post_summary.csv`、`02_doubletfinder_params.csv`、`02_resolution_grid.csv`、`02_celltype_proportion.csv`、`02_cytotrace_by_celltype.csv`、`Fig2*` |
| 03 | `03_QC_lactate_score.csv`、`processed/neo_neoplastic.rds`、`FigS1_*` |
| 04 | `04_QC_scenic.csv`、`Fig3A–H`、`logs/04_0{1,2,3}_*.log` |
| 05 | `05_monocle_state_*_crosstab.csv`、`05_BEAM_branch2_all.csv`、`05_statedependent_GO_BP.csv`、`FigS2*` |
| 06 | `06_rctd_weights_<GSM>.csv`、`06_nhood_enrichment_zscore_<GSM>.csv`、`06_nn_distance_summary.csv`、`06_spot_counts_per_slice.csv`（**11,970 口径的唯一证据，D8**）、`Fig4*` |
| 07 | `07_cellchat_*.csv/rds`、`07_spatial_scheme_AB_comparison.csv`、`07_threshold_sensitivity.csv`、`Fig2C–D`、`Fig4D–E` |
| 08 | `08_Fig3D_NFATc4_survival_summary.csv`、`08_FigS4B/C_*.csv`、`08_qc_cutpoint_*.csv`、`Fig3D`、`FigS4` |
| 09 | `results/ml/*`（`deg_clusterA_vs_B.csv`、`candidate_genes.csv`、`mime1_cindex_matrix.csv`、`auc_timedep.csv`、`dca_1yr_*.csv`、`univariate_cox_meta.csv`、`nested_cv_cindex.csv`）、`Fig7*`、`FigS5` |
| 10 | `Fig8A–E_*`、`Fig8A–E` |
| 11 | `Fig9A–D_*`、`raw/drug_response/**`、`Fig9A–F` |
| 12 | 原始实验数据（qPCR/WB/流式/EdU/Transwell）+ `results/Fig5*/Fig6*` 统计表 |
| 13 | `logs/DECISIONS.md`（决策流水）、gap 表 |
| 14 | 无产物（知识性文档） |

---

## 4. 术语表（本项目唯一口径）

| 术语 | 定义（**唯一口径**） |
|---|---|
| 10 大类细胞标签 | `Neoplastic`、`Macrophage`、`T cell`、`Astrocyte`、`Oligodendrocyte`、`Stem cell`、`Neuron`、`Endothelial`、`ILC`、`Fibroblast`（**大小写/单复数以此为准**；散文可写原文措辞，但代码与图表必须用该 level） |
| 基因 / 蛋白命名 | 人 **`NFATC4`**（别名 `NFAT3`）；蛋白 **NFATc4**；小鼠 **`Nfatc4`**；SPP1 同理（人 `SPP1` / 小鼠 `Spp1`） |
| 高/低乳酸分组 | `lactate_group`：332 LRGs 的 ssGSEA 分数在 **neoplastic 细胞内的中位数二分**（`Low`/`High`）；AUCell/UCell 版本仅作敏感性分析（D14） |
| **NFATc4⁺ 肿瘤细胞** | ① 单细胞层主口径：`NFATC4` regulon AUC 在 neoplastic 内的 **Q75** 以上；② 降级口径（无 regulon 时）：`NFATC4` 表达在 neoplastic 内 Q75 以上；**③ 空间层**：`NFATC4` 表达 ≥ Q75 **且** RCTD 中 neoplastic 权重 ≥ 0.5 |
| **SPP1⁺ 巨噬细胞** | 同构定义：单细胞层为 `SPP1` 表达（或程序分数）在 `Macrophage` 内 Q75 以上；空间层为 `SPP1` 表达 ≥ Q75 **且** macrophage 权重 ≥ 0.5。**定义中不得使用将用于下游结论的同一基因做循环论证**（07 号文档 §4.3 提供正交定义作敏感性分析） |
| 空间区域三分法 | `high-lactate` / `low-lactate` / `non-tumor`：按 spot 的乳酸分数与 **tumor_weight** 划分，阈值写入 `logs/`（文献未报告） |
| 风险分组 | `risk group`：09 号模型风险分数的**中位数二分**（`High`/`Low`）；最优截断点仅作探索性分析 |
| 队列名 | `TCGA`、`CGGA693`、`CGGA325`、`CGGA301`、`GSE108474`、`GSE4412`、`meta`（**写作图上标注：CGGA325 为原文 CGGA625 的替代**） |

---

## 5. 环境与版本钉死

| 组件 | 版本 | 出处 / 说明 |
|---|---|---|
| R | **4.4.2** | 原文指定 |
| Seurat | **5.3.0** | 原文指定；统一使用 5.x API（`FindVariableFeatures`/`DimPlot`/`layer=`） |
| Bioconductor | 3.20 | 与 R 4.4 配套；**完整逐包版本表以 `01_环境搭建与数据获取.md` 为唯一出处** |
| Python | 3.10 | pySCENIC / squidpy / CellTypist / scMetabolism 侧 |
| motif 数据库 | `hg38__refseq-r80__10kb_up_and_down_tss`（feather 与 tbl 同 release） | 版本记录进契约表 |
| 参考基因组 / 注释 | hg38 + GDC 现算（MAF 与 CNV 同版本） | D7 |

**偏离登记**：见 D 表；任何主版本偏离必须写入 `logs/DECISIONS.md` 并在方法学声明。

---

## 6. 统计与随机性（项目特有补充）

| 项 | 本项目取值 |
|---|---|
| 随机种子 | 全局常量 **`20260101`**（R `set.seed(20260101)`；Python `np.random.seed(20260101)`） |
| 乳酸分组阈值 | **中位数**（主）；三分位 / 最优截断仅作敏感性分析 |
| 差异 TF 检验 | 组间 Wilcoxon + BH 校正；报告效应量与显著细胞比例 |
| 空间置换 | **保持空间自相关结构的零模型**（环面平移）；置换次数写入 `logs/` |
| CellChat 置换 | 粗筛 `nboot = 100` → 定稿 `nboot = 1000`（D13） |
| 生存分析 | 主分析连续 Cox + 中位数二分 KM；最优截断点仅探索（并报告 I 类错误膨胀模拟） |
| 多重比较 | 组学层默认 BH(FDR)；原文 "P<0.01" 双口径（raw/adj）复现 |
| 体外实验 | 3 次独立生物学重复，mean ± SD，先 Shapiro-Wilk 再选 t/ANOVA |

其余（检验选择、显著性标注、图表规范）直接遵循 L0 §7。

---

## 7. 决策登记表（D 表）

> **本项目的 D 表权威位置 = `13_复现质控清单、排错与批判性评述.md` §八**（L0 §10.1 允许的"二选一"方案）。
> 本节仅维护**索引与 ★ 标记**；修改 D 表请改 13 号文档，并同步本节。

| # | 决策项 | 一句话 | ★ |
|---|---|---|---|
| D1 | 单细胞对象命名 | 统一 `processed/seurat_gbm_qc_seurat5.rds` | ★ |
| D2 | 目录前缀与 bulk 产物名 | `raw/`、`processed/`；`bulk_<cohort>.rds` | ★ |
| D3 | meta-cohort 构成 | **M1：不含 TCGA**（含则禁称 external validation） | ★ |
| D4 | 模型产物出口名 | `results/ml/stepcox_forward_coefficients.csv` 等 | |
| D5 | "CGGA625" 处理 | 用 **CGGA325** 替代并标注 | ★ |
| D6 | Transwell 板型 | 主用 **Corning 3422**（体积自洽优先） | ★ |
| D7 | 突变/CNV 参考基因组 | hg38 + GDC 现算（MAF 与 CNV 同版本） | |
| D8 | "11,970 spot" 口径 | 逐切片列表为准；11,970 标"原文口径" | ★ |
| D9 | 免疫浸润一致性判据 | ≥2 算法同向（主）/ ≥3（强证据） | |
| D10 | RSF `nodesize` | 主分析 15；敏感性 5/10/15/25 | |
| D11 | `NFATC4` 循环论证检查 | 告警+落盘（不硬中断）+ 剔除后重打分敏感性 | ★ |
| D12 | RCTD `doublet_mode` | 主 `full`，敏感性 `doublet`；图注注明原文未报告 | |
| D13 | CellChat 置换算力 | `nboot=100` 粗筛 → 1000 定稿 | |
| D14 | 乳酸打分双口径 | ssGSEA 主 + AUCell 平行验证，报 Spearman 与分组 Jaccard | ★ |

---

## 8. QC 门槛

- **全局门槛**：`00_目录与复现总览.md` §五（13 条，含 12,618 细胞 / 24 cluster / 10 类 / 6 regulon 模块 / 11,970 spot / C-index 0.87·0.64·0.62 / 5 年 AUC 0.976·0.905 / pooled HR 1.82 / 97.7% 等）。
- **逐项门槛**：`13_复现质控清单、排错与批判性评述.md` §一（32 项 QC 表，含允许波动与 Go/No-Go 判据）。
- **结题**：13 号 §七 gap 表（逐模块填"原文值 / 我方复现值 / 是否一致 / 差异原因 / 处理方式"）。

---

## 9. 偏离声明（对 L0 非红线规则的偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| L0 §3.2 禁空格/中文 | **原文 PDF 与 `NN_*.md` 文档名保留空格与中文** | L0 §3.2 第 6 条已明文豁免（可读性优先） | 无 |
| L0 §10.1 D 表位置 | D 表置于质控文档（13 号 §八）而非本文件 | L0 允许"二选一"，13 号已被大量交叉引用，避免大规模改链 | 需同步本文件 §7 索引 |

其余无偏离。

---

## 10. 修订记录

| 版本 | 日期 | 修订内容 | 说明 |
|---|---|---|---|
| v1.1 | 2026-10-16 | 跟随 L0 升版至 v2.1（新增 §1.4 脚本分层：可复用脚本放仓库 `tools/`，项目专用脚本放本目录 `code/`）；本文件内容无实质变化 | 本版本 |
| v1.0 | 2026-10-16 | 由项目内 `REPRO_SPEC.md`（旧 v1.0，项目级）改写为本 `PROJECT_SPEC.md`：剥离通用条款（上移至仓库根 `../REPRO_SPEC.md` L0 v2.0），保留并新增本项目特有内容（模块编号表、产物契约、**术语表（阳性判定唯一口径）**、D 表索引、偏离声明） | 本版本 |
