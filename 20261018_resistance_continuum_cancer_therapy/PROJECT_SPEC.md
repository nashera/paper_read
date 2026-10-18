# PROJECT_SPEC —— 耐药连续谱（Resistance Continuum）复现项目规范

> **上位规范**：`../REPRO_SPEC.md`（文献复现项目统一规范 v2.1，层级 L0，全仓库强制）
> **本文件层级**：**L1**（L0 在本项目的实例化；不得违背 L0 红线）
> **本项目**：França GS, Baron M, King BR, Bossowski JP, Bjornberg A, Pour M, Rao A, Patel AS, Misirlioglu S, Barkley D, Tang KH, Dolgalev I, Liberman DA, Avital G, Kuperwaser F, Chiodin M, Levine DA, Papagiannakopoulos T, Marusyk A, Lionnet T, Yanai I. *Cellular adaptation to cancer therapy along a resistance continuum.* **Nature (2024) 631:876–890.** DOI: [10.1038/s41586-024-07690-9](https://doi.org/10.1038/s41586-024-07690-9)
> **项目目录**：`20261018_resistance_continuum_cancer_therapy/`
> **启动日期**：2026-10-18　**本文件版本**：**v1.0**（对应 L0 v2.1）
>
> **效力**：L0 `../REPRO_SPEC.md` > **本文件 L1** > 10 号文档 §六 QC 门槛（L2）> 各模块文档（L3）。

---

## 1. 项目定义

### 1.1 一句话目标

在公开的 GEO 处理数据（GSE247691 SuperSeries）与作者代码基础上，**逐图重建**原文证据链：剂量递增产生耐药连续谱 → scRNA-seq 定义 5 个细胞状态与 6 个基因模块 → inferCNV/WES 揭示 CNA 与状态的对应 → ATAC-seq 揭示染色质渐进重塑与转录因子 motif 转换 → PDX 体内复现持久/耐药状态 → 代谢 CRISPR 筛选与代谢组锁定 NRF2/核苷酸代谢依赖 → CB-839 联合用药增敏；并逐条登记文献未写清的参数（D 表）。

### 1.2 原文核心结论（复现时逐条验证）

1. 对 BRCA2 缺陷 HGSOC 细胞系 Kuramochi 做 1→320 μM 奥拉帕利剂量递增，得到 9 个适应群体（T1–T320）；IC50 由对照约 2 μM 升至 T320 约 60 μM（原文 Fig. 1b；属"部分耐药/耐受"而非法完全耐药）。
2. 相同策略在 Ovsaho、COV362（奥拉帕利）、A375（达拉非尼）、PC9（奥希替尼）中重现（Extended Data Fig. 1）。
3. 单细胞转录组沿连续谱出现 5 个主要状态 S1–S5，其频率沿适应线改变；状态可在不同适应线间共享（Fig. 1d–f）。
4. 6 个基因模块 A–F 逐步组装：A（去分化/干扰素）、B（混合上皮–间质）、C（ECM/EMT/缺氧）、D（缺氧/EMT/TGFβ/干性）、E/F（NRF2/ATF4、糖酵解、PPP、核苷酸、铁死亡/铁稳态）（Fig. 2a–c）。
5. EMT/干性程序**使能**适应而非直接等同于耐药：Slug 过表达仅提供部分耐药并加速后续适应（Fig. 2d、Extended Data Fig. 6）。
6. 适应状态与拷贝数改变（CNA）相关：chr7/chr9 区域保守，CDKN2A/B 缺失；未见继发 BRCA1/2 突变；不同细胞系 CNA 谱系特异（Fig. 3a、Extended Data Fig. 7）。
7. ATAC-seq 显示染色质可及性沿连续谱累积变化：谱系因子 SOX17/PAX8/WT1 可及性下降，应激调控因子 AP1/NRF2/ATF4 上升；撤药 120 天的 T320W 部分可逆，再暴露可恢复可及性（Fig. 3b–f、Extended Data Fig. 8）。
8. BRCA2 突变 HGSOC PDX 中，剂量梯度产生持久（persistence）状态，三轮治疗产生稳定耐药（resistance）；耐药肿瘤中出现与体外适应模块 C/D/E 重叠的 cluster 9（Fig. 4、Extended Data Fig. 9）。
9. 代谢 CRISPR 负筛显示沿连续谱逐步获得 NRF2/PPP/核苷酸/铁稳态依赖；代谢组显示 T320 核苷酸及相关前体升高；奥拉帕利 + GLS 抑制剂 CB-839 联合使最适应系 IC50 相对单药下降约 6 倍（Fig. 5、Extended Data Fig. 10）。

### 1.3 交付物清单

| 类别 | 文件 | 状态 |
|---|---|---|
| 原文 | `Cellular adaptation to cancer therapy along a resistance continuum (Nature 2024, 631-876).pdf` | ✅ 已归档 |
| 项目规范 | **`PROJECT_SPEC.md`（本文件，L1）** | ✅ v1.0 |
| 通用规范 | `../REPRO_SPEC.md`（L0 v2.1） | ✅ 已存在 |
| 总览与精读 | `00_目录与复现总览.md` | ✅ |
| 模块文档 | `01`–`11` 共 11 篇 | ✅ |
| 代码 | `code/*.{R,py,sh}`（编号前缀与文档一致） | ⏳ S6 实现 |
| 产物 | `raw/ processed/ results/ figures/ logs/` | ⏳ S6 生成 |

### 1.4 模块编号表

| 编号 | 文档 | 模块 | 复现类型 |
|---|---|---|---|
| 00 | 目录与复现总览 | 文献精读、证据链、数据资产、全局 QC、文档索引 | 索引 |
| 01 | 数据与代码可得性、环境搭建 | GEO/代码库盘点、下载、环境、目录初始化 | 工程 |
| 02 | 剂量递增实验与耐药表型 | 剂量递增、IC50、persister、恒量对照、恒量延长 | 湿实验 SOP |
| 03 | 单细胞转录组 QC 整合与细胞聚类 | inDrop/10X 数据、QC、G1 选择、DrImpute、聚类、DE | 计算 |
| 04 | 细胞状态与基因模块识别 | S1–S5、模块 A–F、GSEA、persister 映射、CD24/CD44 | 计算 |
| 05 | 拷贝数推断与全外显子分析 | inferCNV、WES/Control-FREEC、CNA–状态/模块重叠 | 计算 |
| 06 | ATAC-seq 染色质可及性与转录因子 motif | 比对/peaks/IDR/DESeq2/PCA/motif/T320W 稳定性 | 计算 |
| 07 | PDX 体内模型与单细胞图谱 | 持久/耐药 PDX 设计、scRNA-seq、cluster 映射、IF | 体内 + 计算 |
| 08 | 代谢 CRISPR 筛选与靶向代谢组 | 文库/筛选/MAGeCK-MLE、通路、LC-MS、CB-839 依据 | 计算 + 湿实验 |
| 09 | 功能验证与药物联用 | SOX17/Slug 扰动、活细胞成像、Transwell、WB/qPCR、CB-839 联用、Bliss | 湿实验 SOP |
| 10 | 复现质控清单、排错与批判性评述 | QC 表（L2）、gap 表、D 表权威位置、批判性评述 | 质控 |
| 11 | 背景知识精讲与方法学原理速查 | 生物学/算法原理、公式、术语、文献导读 | 知识 |

---

## 2. 目录与命名（本项目的具体化）

### 2.1 数据源细分

```
raw/
├── geo/
│   ├── GSE206125/      # Kuramochi scRNA-seq（inDrop，C–T320，10 样本）
│   ├── GSE247684/      # A375 + PC9 scRNA-seq（10X CellPlex，6 样本）
│   ├── GSE247685/      # Kuramochi ATAC-seq（bigwig + consensus peaks，28 样本）
│   ├── GSE247686/      # Kuramochi 代谢 CRISPR 筛选（文库 + counts，8 样本）
│   ├── GSE247687/      # Kuramochi WES（CNV + SNV 表，4 样本）
│   ├── GSE247688/      # Kuramochi persisters scRNA-seq（HTO，10 样本）
│   ├── GSE247689/      # Ovsaho + COV362 scRNA-seq（10 样本）
│   └── GSE247690/      # PDX scRNA-seq（persistence + resistance，17 样本）
├── metabolomics/ST003102/   # Metabolomics Workbench 研究 ST003102
└── reference/          # hg38/GRCh38、Ensembl 93、MSigDB、HOMER motifs、blacklist
```

### 2.2 命名特例

| 情形 | 规则 |
|---|---|
| 原文 PDF | 保留原题全称，冒号改连字符（`631:876` → `631-876`），允许空格（L0 §3.2 第 6 条豁免） |
| 适应线标签 | Kuramochi：`C, T1, T2.5, T5, T10, T20, T40, T80, T160, T320`；Ovsaho/COV362：**`C, T1, T2.5, T5, T10, T20, T40`（6 条适应线 + C；scRNA-seq 只测 C、T5–T40）**；A375：`C, T0.03, T0.06, T0.125, T0.25, T0.5, T1, T2, T4`；PC9：`C, T0.008, T0.016, T0.032, T0.075, T0.15, T0.3, T0.6, T1.2` |
| persister | `P10`、`P320`（Kuramochi）；`P1.2`（PC9） |
| 撤药系 | `T320W`（撤药 120 天）；恒量延长系 `T5-E`、`T20-E`、`T80-E` |
| 细胞状态 | `S1`–`S5`（大写 S + 数字） |
| 基因模块 | `A`–`F`（单字母，模块名不得写成 `Module1` 等变体） |
| 出图 | `Fig<N><L>_<语义>.pdf/png` 与 `FigED<N><L>_<语义>.pdf/png`（ED = Extended Data） |
| 模块私有产物 | `<NN>_` 前缀，如 `03_kura_qc_pre.csv`、`06_diff_peaks_T320.csv` |

### 2.3 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `T025`、`T0250`（A375） | `T0.25` | Extended Data Fig. 3i vs Fig. 4 图注不一致 |
| `T05`、`T0500`（A375） | `T0.5` | 同上 |
| `T003`、`T006`、`T0125`（A375） | `T0.03`、`T0.06`、`T0.125` | 同上 |
| `T0008` … `T1.2`（PC9） | `T0.008`、`T0.016`、`T0.032`、`T0.075`、`T0.15`、`T0.3`、`T0.6`、`T1.2` | Extended Data Fig. 3/4 |
| `T0300`、`T0600`（PC9 图注） | `T0.3`、`T0.6` | Extended Data Fig. 4 |
| `P1.2` vs `persister PC9` | `P1.2` | Extended Data Fig. 4 |
| `Module A` … `Module F` | `A` … `F` | 正文 Fig. 2 |
| `resistant line` / `adapted line` | `adapted line`（适应性群体）；`resistant` 仅用于 PDX 三轮治疗后的稳定耐药 | 术语表 §5 |
| `T320-D` / `T320D` | `T320D`（加药状态） | Fig. 3f、Fig. 5 |
| `C-Slug0` / `C-Slug1` | `C-Slug0`（无 dox）/ `C-Slug1`（1 μg/ml dox） | Fig. 2d |
| `CD`（control on drug） | `CD`（C + 1 μM 奥拉帕利） | Fig. 5a |

---

## 3. 产物契约表

> 纪律：**每个产物只有一个生产者**（L0 §4）。改格式/名 = 改本表 + 通知全部消费方 + 记 `logs/DECISIONS.md`。

### 3.1 原始数据（只读）

| 产物 | 生产者 | 格式 / 结构 | 消费方 | 脚本 |
|---|---|---|---|---|
| `raw/geo/GSE206125/GSE206125_kuramochi_all_C_to_T320.csv.gz` | 01（下载） | **基因（symbol）× 细胞** counts；行=A1BG/A2M 等 symbol，列=barcode | 03 | `code/01_download_geo.sh` |
| `raw/geo/GSE206125/GSE206125_kuramochi_C_to_T320_metadata_g1.tsv.gz` | 01 | 列 `cell, nGene, nUMI, orig.ident, percent.mito, State`；`State ∈ {C,S1–S5}`（**GEO 已给出状态标签，可作验证**）；无 phase/ribo | 03/04 | 同上 |
| `raw/geo/GSE247685/*.bw` + `*_consensus_peaks_kuramochi_all.tsv.gz` | 01 | **consensus 文件 = chr/start/end + 20 个样本计数列**（C、T5–T320、P10、P320，各 2 重复）；6 个 bigWig 为系列级轨道 | 06 | 同上 |
| `raw/geo/GSE247686/*sgRNA*.tsv.gz` | 01 | sgRNA 序列库 + per-sample counts | 08 | 同上 |
| `raw/geo/GSE247687/*.tsv.gz` | 01 | WES CNV（segment）与 SNV 表 | 05 | 同上 |
| `raw/geo/GSE247688/`、`GSE247689/`、`GSE247690/` | 01 | counts + metadata（+ HTO barcode / antibody ref） | 03/04/07 | 同上 |
| `raw/metabolomics/ST003102/*` | 01 | 代谢物强度矩阵（若 WB 不可直接下载，则以原文补充表替代，见 D8） | 08 | `code/01_download_metabolomics.sh` |

### 3.2 核心中间对象

| 产物 | 生产者 | 格式 / 结构（下游必须知道的一切） | 消费方 | 脚本 |
|---|---|---|---|---|
| `processed/03_kura_qc.rds` | 03 | `list(mat = 基因×细胞 counts, meta = 细胞元数据)`；`meta` 列含 `cell, nGene, nUMI, orig.ident, percent.mito, State`（GEO 原始）+ 本项目新增 `sample`（= `orig.ident`）、`perc_mito`、`perc_ribo`、`n_counts`、`n_genes`；只保留 QC 通过细胞（GEO 已为 G1 口径）；**phase 在 03 的 Seurat 对象中，不在此 RDS** | 04/05 | `code/03_01_load_qc.R` |
| `processed/03_kura_seurat_list.rds` | 03 | **按样本分别聚类**的 Seurat 对象 list；每个对象含 `seurat_clusters`、PCA、UMAP | 03 出图 | `code/03_04_cluster.R` |
| `results/03_line_clusters.tsv` | 03 | 列 `cell, sample, cluster`（每线独立 Louvain 簇） | 04 | `code/03_04_cluster.R` |
| `processed/03_kura_seurat_g1.rds` | 03 | 合并对象：RNA assay + `line_cluster` 元数据；用于 04/05（不用于跨样本 UMAP） | 04/05 | `code/03_04_cluster.R` |
| `processed/04_state_assignment.tsv` | 04 | 列 `cell, sample, seurat_clusters, state, published_state`（本项目 `state ∈ S1–S5`；`published_state ∈ C,S1–S5` 为 GEO 原始标签，仅作验证） | 05/06/07 映射 | `code/04_01_states.R` |
| `results/04_gene_modules.csv` | 04 | 列 `module, gene, signature`（module ∈ A–F；基因可跨模块） | 05/07/08 | `code/04_02_modules.R` |
| `results/04_state_frequency.csv` | 04 | 列 `sample, state, n_cells, fraction` | 00/10 QC | 同上 |
| `processed/05_infercnv.rds` | 05 | inferCNV 对象（denoised matrix，基因按染色体排序） | 05 出图 | `code/05_01_infercnv.R` |
| `results/05_cna_segments.csv` | 05 | 列 `sample, chr, start, end, log2_ratio, call` | 05/10 | 同上 |
| `results/06_diff_peaks.csv` | 06 | 列 `peak_id, chr, start, end, log2FoldChange, padj, contrast, direction` | 06 motif/T320W | `code/06_01_load_counts.R` |
| `results/06_motif_enrichment.csv` | 06 | 列 `motif, direction, n_peaks, p, p_adj` | 06/10 | `code/06_04_motif.sh` |
| `results/07_pdx_clusters.csv` | 07 | 列 `cell, experiment, sample, cluster` | 07/10 | `code/07_01_pdx.R`（脚本末尾写出） |
| `results/08_crispr_gene_scores.csv` | 08 | 列 `gene, sample, beta, p, fdr`（MAGeCK-MLE） | 08 | `code/08_01_mageck.sh` |
| `results/08_metabolomics_diff.csv` | 08 | 列 `metabolite, comparison, log2FC, p, p_adj` | 08/10 | `code/08_05_metabolites.R` |
| `results/09_ic50_summary.csv` | **02（唯一）** | 列 `line, drug, ic50, sd, n`；单位 μM；6–9 天（卵巢）/4 天（A375）/5 天（PC9）口径分列 | 00/09/10 | `code/02_01_ic50.R` |

### 3.3 图表产物

| 图 | 生产者 | 输入 | 脚本 |
|---|---|---|---|
| `figures/Fig1b_IC50_all_lines.pdf/png` | **02（唯一）** | 09_ic50_summary.csv | `code/02_01_ic50.R` |
| `figures/Fig1d_UMAP_by_line.pdf/png` | 03 | `processed/03_kura_seurat_list.rds` | `code/03_04_cluster.R` |
| `figures/Fig1e_state_dendrogram.pdf/png` | 04 | 04_state_assignment.tsv | `code/04_01_states.R` |
| `figures/Fig1f_state_frequency.pdf/png` | 04 | 04_state_frequency.csv | 同上 |
| `figures/Fig2a_gene_modules_heatmap.pdf/png` | 04 | 04_gene_modules.csv | `code/04_02_modules.R` |
| `figures/Fig2b_module_gsea.pdf/png` | 04 | 04_gene_modules.csv | `code/04_03_gsea.R` |
| `figures/Fig2c_persister_overlap.pdf/png` | 04 | 03/04 产物 | `code/04_04_overlap.R` |
| `figures/Fig3a_cna_heatmap.pdf/png` | 05 | 05_infercnv.rds | `code/05_02_plot_cna.R` |
| `figures/Fig3c_atac_pca.pdf/png` | 06 | 公开 consensus counts（20 样本） | `code/06_02_pca.R` |
| `figures/Fig3e_motif_bars.pdf/png` | 06 | 06_motif_enrichment.csv | `code/06_04_motif.sh` |
| `figures/Fig3f_T320W_heatmap.pdf/png` | 06 | `processed/06_t320w/bigwig` + `results/06_motif/peaks_T320_vs_C.bed` | `code/06_06_t320w.sh` |
| `figures/Fig4d_pdx_umap.pdf/png` | 07 | 07_pdx_clusters.csv | `code/07_01_pdx.R` |
| `figures/Fig4e_pdx_module_overlap.pdf/png` | 07 | 07 + 04 产物 | `code/07_02_overlap.R` |
| `figures/Fig5b_crispr_heatmap.pdf/png` | 08 | 08_crispr_gene_scores.csv | `code/08_03_gene_scores.R` |
| `figures/Fig5e_cb839_ic50.pdf/png` | 09 | `09_viability_combination.csv` / `09_ic50_combination_ratio.csv` | `code/09_05_combination.R` |

---

## 4. 术语表（本项目唯一口径）

| 术语 | 定义（唯一口径） |
|---|---|
| **resistance continuum（耐药连续谱）** | 由一系列细胞状态转换构成、伴随适应度逐步升高的治疗适应轨迹；不是单一开关，也不是预先存在的完全耐药克隆 |
| **adapted line（适应线）** | 剂量递增终点可在此剂量下达到汇合的群体；命名 `T<剂量>`，剂量单位 μM（A375/PC9 为各自的 μM） |
| **persister（持久细胞）** | 单次高剂量（>IC90）短期暴露后仍贴壁存活的细胞；`P10`/`P320`/`P1.2`；表型可逆、无显著 IC50 升高 |
| **resistant（耐药）** | 仅用于 PDX 三轮连续治疗后、再治疗无反应的肿瘤（R1/R2/R3）；体外 adapted line **不**称 resistant |
| **state（状态）** | 由跨样本亚群平均表达谱的 Spearman 相关 + average-linkage 层次聚类（k=5）定义的细胞状态 S1–S5 |
| **module（模块）** | 按 Neftel 方法识别的基因表达程序 A–F；模块间**不要求互斥**（一个基因可属于多个模块） |
| **fitness（适应度）** | 在给定药物剂量下的相对增殖/存活能力；操作化为 IC50 与相对细胞数（CellTiter-Glo） |
| **IC50** | 6–9 天（卵巢癌）或 4 天（A375）/5 天（PC9）处理后，相对 DMSO 对照的 50% 抑制浓度；非线性回归拟合 |
| **G1 细胞主口径** | 所有单细胞分析（状态、模块、CNV）仅使用 CellCycleScoring 判定为 G1 的细胞；全细胞仅作敏感性分析 |
| **CNA call** | inferCNV 的 denoised 表达相对参考（未处理 C）的偏离；`gain`/`loss` 阈值以 inferCNV 默认 HMM 为准 |
| **高/低可及性 peak** | 相对 C，DESeq2 `p_adj < 0.001` 且 `|log2FC| ≥ log2(3)` |
| **T320W** | T320 撤药 120 天的群体；用于检验染色质状态的稳定性与"记忆" |
| **dose-escalation** | 每次达到 >70% 汇合后传代并以倍增剂量继续的适应方案 |
| **constant-dose** | 在固定剂量下持续培养相同累计时间的对照方案（`T5-E`/`T20-E`/`T80-E`） |
| **随机种子** | 全项目 `20261018` |

---

## 5. 环境与版本钉死

| 组件 | 版本 | 出处 / 说明 |
|---|---|---|
| R | 4.4.2 | 项目环境；原文未报告 R 版本 |
| Seurat | **4.4.0** | 原文使用 4.0.5 与 2.3.2；选 4.x 末版以最大兼容旧 API（D1） |
| Bioconductor | 3.20 | 与 R 4.4 配套 |
| inferCNV | 1.4.0 | 原文指定 |
| ComplexHeatmap | 2.11.1 | 原文指定（现代版本可用，图样式可能微差） |
| fgsea | 1.16 | 原文指定 |
| GeneOverlap | 1.26 | 原文指定 |
| UCell | 2.0.0 | 原文指定 |
| regioneR | 1.10.0 | 原文指定 |
| Mfuzz | 2.5.0 | 原文指定 |
| MAST | 1.18+ | Seurat `test.use="MAST"` |
| Python | 3.10 | MAGeCK、deepTools、统计脚本 |
| MAGeCK | 0.5.8 | 原文指定 |
| Bowtie2 | 2.3.4.1 | 原文指定 |
| MACS2 | 2.2.4 | 原文指定 |
| HOMER | 4.15 | 原文指定 |
| deepTools | 3.1 | 原文指定 |
| DESeq2 | 1.18.1 | 原文指定；现版本结果趋势一致 |
| IDR | 2.0.3 | 原文指定 |
| TrimGalore | 0.5.0 | 原文指定 |
| Picard | 2.18.2 | 原文指定 |
| SAMtools | 1.9 | 原文指定 |
| Control-FREEC | 11.x | 原文未给版本 |
| CellRanger | 3.1 / 7.0.0 | 原文指定（仅重跑 raw fastq 时需要） |
| 参考基因组 | hg38/GRCh38；Ensembl 93 | 原文指定 |
| motif 库 | HOMER known motifs（hg38） | 原文指定工具；具体库版本未报告 |

**环境偏离登记**：任何主版本偏离写入 D 表与 `logs/DECISIONS.md`；Seurat 4.4.0 与原文 4.0.5 的差异见 D1。

---

## 6. 统计与随机性（项目特有补充）

| 项 | 本项目取值 |
|---|---|
| 随机种子 | **`20261018`**（R `set.seed(20261018)`；Python `np.random.seed(20261018)`；MAGeCK 无随机性） |
| 差异表达 | MAST；`logfc.threshold=0.25`、`min.pct=0.2`、**`min.cells.feature=10`**（Seurat 4.4 参数名；原文写作 `min.cells.genes`）；BH 校正 |
| 状态定义 | Spearman 相关 + average linkage；k=5（敏感性 k=4/6） |
| 模块识别 | **`avg_log2FC > log2(1.8)≈0.85`**（原文"1.8 average fold change"）；`p_adj < 0.001`；>50 基因；Jaccard >0.75 合并；≥3 signatures |
| 模块 GSEA | fgsea；10,000 次置换；P<0.05 且 FDR<0.25（Fig. 2b 用 P<0.01、FDR<0.25） |
| CNA 重叠 | regioneR；10,000 次置换；双侧 |
| ATAC 差异 peak | DESeq2；`p_adj<0.001` 且 FC≥3 |
| ATAC motif | HOMER；one-sided hypergeometric；BH；报告 P<1e-5 的 TF |
| ATAC–CNA 重叠 | regioneR；1,000 次置换；双侧 |
| CRISPR 筛选 | MAGeCK-MLE；负 β、P<0.05、FDR<0.25；score diff >1 SD |
| 代谢组 | 双侧 t 检验；P<0.05；S/N≥3 且下限 10,000 |
| 细胞实验 | 3 次独立生物学重复；mean±SD；先 Shapiro-Wilk，再 t/ANOVA；IC50 用非线性回归 |
| 多重比较 | 组学层默认 BH(FDR)；原文只给 P 值处同时报告 raw 与 adj |
| 图表 | 每图可追溯"脚本 + 输入产物 + 参数"；图注写 n、检验、p 形式 |

---

## 7. 决策登记表（D 表索引）

> **权威位置**：`10_复现质控清单、排错与批判性评述.md` §八。本节只维护索引；★ = 改变结论口径。

| # | 决策项 | 一句话 | ★ |
|---|---|---|---|
| D1 | Seurat 版本 | 主口径 Seurat 4.4.0；Seurat 5 仅作敏感性 | ★ |
| D2 | 数据来源 | 优先使用 GEO 处理后 counts/metadata；raw fastq 重跑为备选 | ★ |
| D3 | 细胞周期 | GEO 已按 G1 过滤；主口径保留全部并另算 phase 复核 | ★ |
| D4 | 状态数 k | 固定 k=5；k=4/6 作敏感性 | ★ |
| D5 | 模块识别实现 | 严格复刻 Neftel 流程；参数见 04 文档；结果不一致时报告而非调参 | ★ |
| D6 | inferCNV 参考 | 未处理 C 为参考；cutoff 0.1；HMM 默认 | |
| D7 | ATAC 差异分析输入 | **公开 consensus 文件已含 20 样本计数**（C、T5–T320、P10、P320）：Fig. 3d/e 与 PCA 直接复现；**Fig. 3f 的 T320/T320W 加药 8 样本未公开**，需 SRA 重跑；6 个 bigWig 仅作信号可视化 | ★ |
| D8 | 代谢组数据 | 优先 ST003102；不可得时用补充表做定性核对，不做定量复现 | ★ |
| D9 | PDX 批次校正 | fastMNN，与原文一致；不额外去批次 | |
| D10 | CRISPR 参考比较 | CD/T10/T320 vs C；T320D vs T320 与 CD；**T40 列存在但非主分析**（敏感性用） | ★ |
| D11 | 体外"耐药"措辞 | adapted line 不称 resistant；仅 PDX 三轮后称 resistant | ★ |
| D12 | A375/PC9 线名 | 统一小数命名（见 §2.3），图注保留原文别名对照 | |
| D13 | 湿实验范围 | 仅写 SOP 与 QC，不执行；计算复现以公开数据为准 | ★ |
| D14 | 代码库使用 | 作者代码仅作参考；缺失脚本按文档重建，不假设未公开中间对象 | ★ |
| D15 | 随机种子 | 20261018 | |
| D16 | 模块数 | 以 Neftel 流程的签名聚类为准；目标 6 个（A–F），允许 5–7 并报告 | ★ |
| D17 | GSEA 阈值 | 主口径 P<0.05、FDR<0.25；Fig. 2b 用 P<0.01、FDR<0.25 | |
| D18 | 多重校正 | 代谢组同时报告 raw P 与 BH；其余默认 BH | |
| D19 | IC50 拟合 | 4PL（`drc::LL.4`，ED50 = 参数 `e`） | |
| D20 | 图注 n | 每图写 n、检验方法、p 值形式 | |
| D21 | 归一化口径 | GEO 给 counts；主口径 Seurat `LogNormalize`（CP10k），**不等于**原文 TPM | ★ |
| D22 | GEO `State` 标签 | 仅作验证/对照（ARI 与频率），不直接复制为本项目结果 | ★ |

---

## 8. QC 门槛

- **全局门槛**：`00_目录与复现总览.md` §六（逐图 Go/No-Go，含原文报告值）。
- **逐项门槛**：`10_复现质控清单、排错与批判性评述.md` §六（模块级 QC 表，含允许波动）。
- **结题**：`10_复现质控清单、排错与批判性评述.md` §十 gap 表。

---

## 9. 偏离声明（对 L0 非红线规则的偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| L0 §1.2 文件名禁空格 | 原文 PDF 保留空格 | L0 §3.2 第 6 条明文豁免（可读性优先） | 无 |
| L0 §1.2 | 项目 slug 含 `resistance_continuum_cancer_therapy`（4 词） | 在 3–6 词范围内 | 无 |
| L0 §10.1 | D 表置于 10 号质控文档而非本文件 | L0 允许"二选一"；便于交叉引用 | 需同步 §7 索引 |

---

## 10. 修订记录

| 版本 | 日期 | 修改内容 | 修改人 |
|---|---|---|---|
| v1.0 | 2026-10-18 | 首版：项目立项、模块划分、产物契约、术语表、D 表索引、环境与 QC 入口 | — |
