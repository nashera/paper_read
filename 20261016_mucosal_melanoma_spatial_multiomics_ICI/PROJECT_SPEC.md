# PROJECT_SPEC —— MuM 空间多组学 ICI 响应 复现项目规范

> **上位规范**：`../REPRO_SPEC.md`（文献复现项目**统一规范（通用）v2.1**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（通用规范在本项目的实例化；不得违背 L0 红线）
> **本项目**：Wang J, Nagarajan P, et al. *Integration of spatial single-cell proteomics and spatial metabolomics reveals tumor microenvironment predictive of immunotherapy response in mucosal melanoma.* **bioRxiv 预印本（未经同行评审）**，2026-04-27. DOI: 10.64898/2026.04.23.720453
> **项目目录**：`/mnt/e/paper_read/20261016_mucosal_melanoma_spatial_multiomics_ICI/`
> **启动日期**：2026-10-16　**本文件版本**：**v1.0**（对应 L0 v2.1）
>
> **效力**：本项目内，L0 `../REPRO_SPEC.md` > **本文件 L1** > 12 号文档 §一 QC 门槛（L2）> 各模块文档（L3）。
> **D 表权威位置**：`12_复现质控清单、排错与批判性评述.md` §七。

---

## 1. 项目定义

### 1.1 一句话目标

在**原始数据不公开**的现实约束下，分三层复现本文：**A 类**（算法与代码，用合成数据验证正确性）→ **B 类**（在公开替代数据上跑通全流程）→ **C 类**（患者级结论，明确标注不可复现并给出数据索取路径）。

### 1.2 核心特殊性：三层复现性框架（本项目一切取舍的依据）

| 层 | 内容 | 可做到 |
|---|---|---|
| **A 算法层** | CN（70 μm/kNN 窗口 + K-means）、空间上下文图（>85% 最少组合）、重心坐标界面、patch 邻近（HDBSCAN+凹包+3NN+5–25 μm 环带）、Delaunay+置换相互作用、PLS-DA、KEGG/MSEA、Spearman+BH、混合效应 | **完全复现**，并以合成数据验证（ARI>0.7、功效>0.8、FPR<0.1、环带误差<15%） |
| **B 流程层** | 上述流水线在公开替代数据（HTAN 多重成像、METASPACE MALDI-MSI）上端到端跑通 | 流程可复现；**结果数字不可与原文比较**（面板/疾病不同） |
| **C 结论层** | 26 例、97 cores 的患者级空间与代谢结论 | **不可复现**；只能索取数据或做同构假说检验 |

### 1.3 模块编号表

| 编号 | 文档 | 模块 |
|---|---|---|
| 00 | 目录与复现总览 | 三层框架、索引、流程图、图表映射、全局门槛 |
| 01 | 数据可得性评估与替代数据策略 | 数据审计、索取话术、替代数据集核实、合成数据生成器、伦理 |
| 02 | 实验平台 SOP（TMA-COMET-MALDI） | TMA、COMET 40-plex seqIF、黑色素漂白、1,5-DAN、timsTOF 参数 |
| 03 | 影像预处理配准与细胞分割 | flat-field、自发荧光、跨模态配准、ROI、核分割与扩张 |
| 04 | 细胞表型定义与注释 | arcsinh、批次校正、8 大类/25 状态、PLS-DA 边界、Voronoi |
| 05 | 细胞邻域 CN 分析与 K-means | 窗口口径、组成向量、k=15、稳定性、SPACEc 替代 |
| 06 | 空间上下文图与邻域界面分析 | >85% 最少组合覆盖、barycentric 投影（≥80%） |
| 07 | Patch 邻近分析 | HDBSCAN、凹包、3-NN、5–25 μm 环带（含密度归一化） |
| 08 | 细胞间相互作用推断与置换检验 | Delaunay、265 px、100/1000 置换、Mann–Whitney、logFC |
| 09 | 临床关联与无监督空间分型 | R/NR 与 INMC、伪重复处理、5 簇、功效分析 |
| 10 | 空间代谢组学分析 | MALDI 处理链、632/82、PLS-DA、MSEA、tryptophan score |
| 11 | 蛋白-代谢物空间耦合与多组学整合 | 配准误差传播、Spearman+BH、空间零模型、Simpson 悖论、Fig. 6 |
| 12 | 复现质控清单、排错与批判性评述 | QC 表、排错、**原文 25 处疏漏**、**D 表**、批判性评述 |
| 13 | 背景知识精讲与方法学原理速查 | 生物学、平台原理、方法学速查、50 篇文献导读 |

---

## 2. 目录与命名（本项目具体化）

### 2.1 目录

```
raw/{comet,maldi,clinical,tma,public_surrogate}/   # 原始与替代数据（只读）
processed/     # 中间对象（parquet/h5ad/npy/zarr/tif）
results/       # 表格结果（csv/tsv/json）
figures/       # 出图（pdf + png，dpi=300）
code/          # 项目编排代码（仅项目特有参数与流程；通用算法调用 tools/）
logs/          # sessionInfo / DECISIONS.md / 校验和
```

### 2.2 命名特例

| 情形 | 规则 |
|---|---|
| 原文 PDF | 保留原题全称（含空格），冒号改连字符 |
| 图表 | 与原文图号对齐：`Fig<N><L>_<语义>.pdf/png`（如 `Fig2D_preR_vs_preNR.pdf`、`Fig4I_hierarchical_clusters.pdf`） |
| 核心/切片占位 | `<core_id>`、`<tma>`、`<GSM>`、`<pair>`（**统一尖括号，禁止 `{}`**） |
| 模块私有产物 | `<NN>_` 前缀（`03_qc_segmentation.csv`、`05_cells_cn.parquet`） |
| 可复用算法 | **放仓库 `../tools/`**（`spatial_cn_kmeans.py`、`patch_proximity.py`、`spatial_interaction_test.py`），项目内不得复制 |

### 2.3 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `processed/03_ilastik_roi_{nickname}.tif` | `processed/03_ilastik_roi_<core_id>.tif` | L1（占位符统一为 `<>`） |
| `processed/05_cn_composition.npz` 与 `06_cn_composition_knn70.npz` 并存 | 前者为 **70 μm 半径版**（主口径 D4），后者为 **70 kNN 版**（敏感性）；文件名须显式带口径标签 | D4 |
| `figures/Fig2D_cn_boxplot.pdf` | `figures/Fig2D_preR_vs_preNR.pdf` + `figures/Fig2D_preNR_vs_postNR.pdf`（两个对比分开出图） | L1 |
| `Fig3A_context_subset_graph.pdf` | 作为 `Fig3A_spatial_context_map.pdf` 的补充件保留，命名前缀 `Fig3A_` 一致即可 | L1 |

---

## 3. 产物契约表

> 纪律（L0 §4）：**每个产物唯一生产者**；改格式/名 = 改本表 + 通知消费方 + 记 `logs/DECISIONS.md`。

### 3.1 核心共享对象（无前缀）

| 产物 | 生产者 | 格式 / 结构（下游必须知道的） | 消费方 |
|---|---|---|---|
| `processed/03_cells_roi.h5ad` | 03 | 细胞表：`core_id, x, y, area, roi_label(tumor/TME), tma`；**坐标为像素，并记录 `pixel_size_um=0.23`** | 04–11 |
| `processed/04_cells_phenotyped.parquet` | 04 | 加 `marker_*`（arcsinh 后）、`celltype`（8 大类）、`phenotype`（25 状态） | 05–09 |
| `processed/05_cells_cn.parquet` | 05 | 加 `cn`（1–15）、`cn_radius`/`cn_knn` 两套标签（口径 D4） | 06–11 |
| `processed/05_cn_composition.parquet` | 05 | CN × 细胞类型组成矩阵（**L1 归一化**） | 06/09 |
| `processed/06_spatial_context.parquet` | 06 | 每细胞的 spatial context（>85% 最少 CN 组合）与层级 | 09/11 |
| `processed/07_patch_labels.parquet` | 07 | 细胞级 patch 标签（-1 = 非 patch） | 08（相互作用分层） |
| `processed/08_edges.parquet` | 08 | Delaunay 边：`cell_i, cell_j, dist_px, dist_um, core_id` | 08–11 |
| `processed/10_metabolites_<core>.parquet` | 10 | MALDI 特征矩阵（RMS 归一化后），索引为像素坐标 | 10/11 |
| `processed/11_cell_metab_link.parquet` | 11 | 细胞/区域 ↔ 代谢物像素的空间连接表 | 11 |
| `results/01_compactsim_*` | 01 | **合成数据（ground truth）**：细胞表 + 像素表 + `01_compactsim_ground_truth.json` | 05–11 的验证基准 |

### 3.2 模块产物索引

> ⚠️ **本表与各模块文档 §⑧「图表对应表」互为校验**：以**模块文档 §⑧ 为全量权威**，本表只列**跨模块被消费**的关键产物。
> 命名不一致时以模块文档为准，并须回改本表（L0 §11.1 变更流程）。**表中名称已于 2026-10-16 按各模块文档实测核对一遍。**

| 模块 | 跨模块消费的关键产物 |
|---|---|
| 01 | `results/01_public_surrogate_candidates.csv`、`01_surrogacy_scores.csv`、`01_metaspace_datasets.csv`、`01_aggregation_bias.csv`、`01_delaunay_unit_of_analysis.csv`、`01_pathway_prior_GSE270464.csv`、`results/01_compactsim_ground_truth.json` |
| 02 | `results/02_melanin_bleach_qc.csv`、`02_dan_vs_analyte_collisions.csv`、`02_fiducial_pairs.csv`、`02_maldi_tic_qc.csv`、`02_registration_tre.csv`、`figures/02_tre_simulation.pdf` |
| 03 | `processed/03_cells_roi.h5ad`、`processed/03_cells_<tma>.parquet`、`processed/03_core_manifest.csv`、`03_roi_mask_<core_id>.tif`、`results/03_qc_segmentation.csv`、`03_registration_sensitivity.csv`、`03_seg_agreement.csv` |
| 04 | `processed/04_cells_phenotyped.parquet`、`results/04_phenotype_fractions.csv`、`04_cluster_vs_manual_ari.csv`、`04_cofactor_sensitivity.csv`、`04_plsda_cv.csv`、`04_suppfig1_*` |
| 05 | `processed/05_cells_cn.parquet`（+ h5ad 镜像，**同由 05 唯一生产**）、`processed/05_cn_composition.parquet`、`results/05_cn_composition_profile.csv`、`05_cn_naming.csv`、`05_cn_frequencies_per_core.csv`、`05_k_sweep_stability.csv` |
| 06 | `processed/06_spatial_context.parquet`、`results/06_context_combination_frequencies.csv`、`06_barycentric_zones_per_core.csv`、`06_barycentric_tau_sensitivity.csv`、`06_cn10_vs_immune_correlation.csv`、`06_cn10_group_assign.csv` |
| 07 | `processed/07_patch_labels.parquet`、`07_patch_anchors.parquet`、`07_patches.parquet`、`07_patch_hull_polygons.parquet`、`07_spanning_points.parquet`、`results/07_proximity_counts_by_radius.csv`、`07_patch_summary.csv`、`07_gradient_test.csv`、`07_hdbscan_param_sensitivity.csv` |
| 08 | `processed/08_edges.parquet`、`08_obs_pairs.parquet`、`08_null_B1000.parquet`、`results/08_interaction_tests.csv`、`08_edges_summary.csv`、`08_core_level_log2FC.csv`、`figures/08_synthetic_fpr_power.pdf` |
| 09 | `processed/09_cn_freq.h5ad`、`results/09_cluster_assignments.csv`、`09_cluster_enrichment.csv`、`09_cluster_ari.csv`、`09_cluster_counts_reconstructed.csv`、`09_cluster_stability.csv`、`09_cn_freq_vs_group.csv`、`figures/09_power_curve.pdf` |
| 10 | `processed/10_maldi_peakmatrix.npz`、`10_maldi_peaks.h5ad`、`10_maldi_pixels_raw.npz`、`10_maldi_core_roi.npz`、`results/10_diff_metabolites_tumor_rms.csv`、`10_diff_metabolites_tme_rms.csv`、`10_key_metabolite_calc_mz.csv`、`10_plsda_q2_perm.csv`、`10_msea_kegg.csv`、`10_tryptophan_score.csv` |
| 11 | `processed/11_linked_grid.npz`、`results/11_corr_partial_{...}.csv`、`11_multiscale.csv`、`11_n_eff.csv`、`11_null_torus.csv`、`11_covariates_protein.csv`、`11_registration_tre.csv` |
| 12 | `logs/DECISIONS.md`、gap 表 |

---

## 4. 术语与判定口径（唯一化）

| 术语 | 定义（**唯一口径**） |
|---|---|
| MuM | mucosal melanoma，黏膜黑色素瘤 |
| R / NR | responder / non-responder；**判定标准按 D10**（原文未报告 → RECIST 1.1 或 iRECIST；post-treatment 用 INMC 病理反应） |
| 三组 | **pre-treatment R、pre-treatment NR、post-treatment NR**（**post-R 缺失**，任何"治疗前后"结论必须声明此限制） |
| CN | cellular neighborhood；**主口径 = 70 μm 半径窗口 + K-means(k=15)**（D4）；kNN70 为敏感性分析 |
| 阳性标签（Ki67⁺/CD163⁺/CD11c⁺/MIF⁺/GZMB⁺/FOXP3⁺/IRF4⁺/Arg1⁺） | **一律按 04 号文档的 marker 规则表**（阈值必须在 `results/04_threshold_robustness.csv` 中记录；**禁止在不同文档中另立阈值**） |
| 像素↔μm | COMET **0.23 μm/px**（主）、MALDI **20 μm/px**；**任何长度/面积参数必须注明单位并断言**（L0 §附录 C） |
| 相互作用阈值 | **265 px ≈ 60.95 μm**（D6），两口径都报 |
| logFC（相互作用） | `log2(观测距离均值 / 置换距离均值)`；**负值 = 观测更近 = 空间富集**（全项目统一） |
| tumor / TME ROI | SOX10 驱动的 ML 分类 → tumor ROI；**其余组织为 TME ROI**（含基质、免疫、血管、正常黏膜） |
| tryptophan score | 三种实现取一致性（D13）：平均 z-score / ssGSEA / 加权和 |
| MSI 置信度 | 本文代谢物注释为 **level 2/3（putative）**；引用具体代谢物时必须带此限定 |

---

## 5. 环境与版本钉死

| 组件 | 版本 | 说明 |
|---|---|---|
| Python | 3.10+ | 原文"所有统计在 Python 中完成"，本项目主语言 |
| 必需包 | numpy, pandas, scipy, scikit-learn, matplotlib, seaborn, h5py, pyarrow, tifffile, zarr | 见 `../tools/requirements.txt` |
| 可选包 | hdbscan, shapely, anndata, squidpy, scimap, cellpose, stardist, pyopenms, malab, statsmodels, pypls | 缺失时按各文档的回退路径 |
| 商业软件（原文） | Visiopharm VIS 2024.07.1.16912、SCiLS Lab 2025b、MetaboScape、FlexImaging 7.0 | **闭源**；开源替代见 D3 |
| 参考库版本 | KEGG release（记录日期）、HMDB 版本、MetaboScape 库版本 | **必须写入 `results/` 元数据** |

---

## 6. 统计与随机性

| 项 | 取值 |
|---|---|
| 随机种子 | `np.random.seed(20260101)`（K-means、HDBSCAN、PLS-DA、置换、聚类、抽样） |
| 组间比较 | Mann–Whitney U（连续/比例，与原文一致）；**主分析用患者级聚合或混合效应（β-二项/GEE）**（D9） |
| 相关 | Spearman + BH（原文口径）；**另加空间零模型**（环面平移）作稳健性 |
| 多重校正 | BH(FDR)；代谢物差异**双口径**（原文 p<0.05&FC>1 + BH 校正版，D11） |
| 置换 | 相互作用：复现 100 次 + 定稿 1000 次（D8）；**必须报告 p 的分辨率下限 1/n_perm** |
| 分类器验证 | 嵌套交叉验证；PLS-DA 报 Q² / 成分数 / 置换检验（D12） |
| 阈值 | p 或 FDR < 0.05 |

其余遵循 L0 §7。

---

## 7. 决策登记表（D 表）索引

> **权威位置 = `12_复现质控清单、排错与批判性评述.md` §七**（**D1–D24**）。本表仅列 ★（改变结论口径）项；新增 D16–D24（265 px 口径、Fig. 4I 的 p 值口径、凹包参数、患者级模型、环带 vs 累积、FOV 定义、细胞尺度常数、CN 标签层级、外部资源计数时效）见该处：

| # | 决策项 | 一句话 | ★ |
|---|---|---|---|
| D1 | 数据获取策略 | 索取与替代并行，不等待 | ★ |
| D2 | 复现目标层级 | 以 A 类为主交付，C 类只做"不可复现声明" | ★ |
| D3 | 商业软件替代 | 开源栈替代并记录偏差 | ★ |
| D4 | CN 窗口口径 | 主用 70 μm 半径，kNN70 作敏感性 | ★ |
| D6 | 265 px 换算 | 0.23 μm/px → ≈60.95 μm，两口径都报 | ★ |
| D7 | 相互作用零模型 | 位置置换为主 + 环面平移/标签置换 | ★ |
| D9 | 统计层级 | 患者级/混合效应为主，core 级仅描述 | ★ |
| D10 | R/NR 判定 | RECIST 1.1/iRECIST + INMC | ★ |
| D11 | 代谢物差异阈值 | 双口径（原文 + BH） | ★ |
| D12 | PLS-DA 验证 | 嵌套 CV + 置换 + Q²/VIP | ★ |
| D13 | tryptophan score | 三实现取一致性 | ★ |
| D15 | 可复用脚本归属 | 三个空间算法放 `tools/` | ★ |

---

## 8. QC 门槛

- **全局门槛**：`00_目录与复现总览.md` §五（8 条原文数值 + **A 类验证门槛**）。
- **逐项门槛**：`12_...md` §一（**40 项** QC 表，含层归属与判据）。
- **A 类验收（本项目的核心判据）**：`tools/spatial_cn_kmeans.py --selftest`（ARI > 0.7）、`tools/spatial_interaction_test.py --selftest`（功效 > 0.8、FPR < 0.1）、`tools/patch_proximity.py --selftest`（环带误差 < 15%）。
- **结题**：`12_...md` §六 gap 表。

---

## 9. 偏离声明（对 L0 非红线规则的偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| L0 §3.2 禁空格/中文 | 原文 PDF 与 `NN_*.md` 保留空格与中文 | L0 §3.2 明文豁免 | 无 |
| L0 §10.1 D 表位置 | D 表置于 12 号文档而非本文件 | L0 允许"二选一"；12 号已被大量交叉引用 | 需同步本文件 §7 索引 |
| L0 §1.4 | 本项目**不新增** `tools/` 脚本，直接调用已有 3 个空间算法脚本 | 复用优先（L0 §1.4 规则 1） | 无 |

---

## 10. 修订记录

| 版本 | 日期 | 修订内容 | 说明 |
|---|---|---|---|
| v1.0 | 2026-10-16 | 首版：三层复现性框架、模块编号表、产物契约表、**术语与判定口径唯一化**、环境与统计取值、D 表索引、A 类验收门槛、偏离声明 | 本版本 |
