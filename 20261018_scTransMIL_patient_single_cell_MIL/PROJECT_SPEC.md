# PROJECT_SPEC —— scTransMIL 患者级-单细胞 MIL 复现项目规范（L1）

> **上位规范**：`../REPRO_SPEC.md`（通用规范 **v2.1**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（L0 在本项目的实例化；不得违背 L0 红线）
> **本项目**：Tang Z, Wang F, Yang F, Song J, Li Y, Zhou J, Song Y, Chen S, Zhu J, You L, Chen CY-C, Yao J. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* **Nature Communications (2026) 17:10800.** DOI: [10.1038/s41467-026-77538-5](https://doi.org/10.1038/s41467-026-77538-5)
> **项目目录**：`20261018_scTransMIL_patient_single_cell_MIL/`
> **启动日期**：2026-10-18　**本文件版本**：**v1.0**（对应 L0 v2.1）
>
> **效力**：L0 `../REPRO_SPEC.md` > **本文件 L1** > 13 号文档 §五 QC 门槛（L2）> 各模块文档（L3）。
> **D 表权威位置**：`13_复现质控清单、排错与批判性评述.md` §八；本文件 §8 只维护索引。

---

## 1. 项目定义

### 1.1 一句话目标

以作者公开代码 `TencentAILabHealthcare/scTransMIL`（Zenodo `10.5281/zenodo.21151067`）为基准，**重训并核验**"rank 编码 + 6 层 Transformer + 门控注意力 + 多实例学习（MIL）"这条"样本级标签 → 单细胞级疾病相关性"的流水线；在泛癌图谱、甲状腺癌、外部 BRCA、AML 四个数据域上复现样本级诊断、单细胞级肿瘤评分、癌种分类、注意力基因发现与 AML 亚型/生存分析，并逐条登记原文未写清、正文与代码不一致的参数（D1–D23）。

### 1.2 原文核心结论（复现时逐条验证）

1. **框架**：把每个样本（或 ≤256 细胞的 sub-sample）当作 MIL 的 **bag**，每个细胞当作 **instance**；只用 bag 级标签训练，即可同时输出样本级诊断、单细胞级 tumor score 与基因级 attention（Fig. 1）。
2. **编码**：每个细胞按表达量做 **rank value encoding**（去掉零表达基因、按表达量排序、映射到 Geneformer 词表），送入 6 层 Transformer，得到 256 维细胞嵌入；rank 编码对绝对表达值不敏感，可缓解批次效应。
3. **聚合**：门控编码器 `h_gate = tanh(W1h) ⊙ sigmoid(W2h)`，tumor score `a = wᵀh_gate`；bag 聚合为 `Σ a_i h_gate_i`（正文 Eq.3；代码实现为 softmax attention 加权，见 D3）。
4. **损失**：bag 级 BCE/CE + **adaptive clustering loss**（按 tumor score 排序取 top-B/bottom-B 作伪标签，B=8，smooth SVM，β=0.7）。
5. **泛癌**：4 M cells / 1,563 samples / 999 donors / 30 cancer types / 104 datasets / 30,000 genes；按 donor 80/10/10 划分，n=7 次重复；scTransMIL 在 Acc/AUC/F1 上优于 KIDA、DEGAS、scIDST、scVI+MIL，推理 245.16 s/bag（DEGAS 720 s）（Fig. 2）。
6. **单细胞**：以 CNV 注释为参照，泛癌单细胞 AUC = **0.71**（scVI+MIL 0.67、scIDST 0.58、DEGAS 0.56）（Fig. 2f）。
7. **癌种**：30 类原发/转移癌种分类；外部 BRCA（37 samples / 40 万细胞）泛化；scIB 批次校正/生物保守评分优于变体（Fig. 3）。
8. **甲状腺早筛**：GSE184362 中 Patient 5 的 paratumor 样本被 scTransMIL 判为阳性（score **0.64**），高 tumor score 集中在癌前上皮细胞（Fig. 2d–e）。
9. **ALL 基因发现**：tumor B 细胞的 top-100 attention 基因包含 CD79A、PIM1，并找回常规 DEA 漏掉的 **BLK**；top-20 通路 FDR<0.001（Fig. 4）。
10. **AML**：GSE116256（35 tumor + 6 normal、30,000 cells、21 donors）；仅用 20 tumor + 3 normal 样本级标签微调，测试 15 + 3；识别 6 个 AML 亚型（HSC-like … Mono-like），HSC-like 生存更差（bulk 671 例，dbGaP）（Fig. 5）。

### 1.3 交付物清单

| 类别 | 文件 | 状态 |
|---|---|---|
| 原文 | `2026_scTransMIL_bridges_patient-level_disease_states_and_single-cell_transcriptomics_for_cancer_screening_and_heterogeneity_inference.pdf` | ✅ 已归档 |
| 项目规范 | **`PROJECT_SPEC.md`（本文件，L1）** | ✅ v1.0 |
| 通用规范 | `../REPRO_SPEC.md`（L0 v2.1） | ✅ |
| 总览与目录 | `00_目录与复现总览.md` | ✅ |
| 模块文档 | `01`–`12` 共 12 篇 | ✅ |
| 质控与批判 | `13_复现质控清单、排错与批判性评述.md`（NN-1） | ✅ |
| 背景速查 | `14_背景知识精讲与方法学原理速查.md`（NN） | ✅ |
| 代码 | `code/*.py/.sh`（编号与文档一致） | ⏳ S6 实现 |
| 运行日志 | `logs/DECISIONS.md`、`logs/13_review_record.md`、`logs/<NN>_pip_freeze.txt` | ⏳ S6 生成 |

### 1.4 模块划分与编号表

| 编号 | 文档标题 | 对应原文位置 | 状态 |
|---|---|---|---|
| 00 | 目录与复现总览（`00_目录与复现总览.md`） | Fig. 1；全文 | ✅ |
| 01 | 环境搭建与数据获取（`01_环境搭建与数据获取.md`） | Data/Code availability；Methods「Data acquisition」 | ✅ |
| 02 | 数据预处理与 MIL 样本构造（`02_数据预处理与MIL样本构造.md`） | Methods「Data preprocessing」；`1-process_v3.py` | ✅ |
| 03 | 模型架构：rank 编码、Transformer、门控与聚合（`03_模型架构_rank编码_Transformer_门控与聚合.md`） | Fig. 1c–d；Methods「Instance/Gated encoder」「Bag aggregator」 | ✅ |
| 04 | 损失函数与训练策略（`04_损失函数与训练策略.md`） | Fig. 1e；Methods「Loss」「Training」 | ✅ |
| 05 | 推理、评估与数据划分（`05_推理_评估与数据划分.md`） | Methods「Inference」「Statistics and reproducibility」 | ✅ |
| 06 | 基线方法与变体复现（`06_基线方法与变体复现.md`） | Methods「Comparison methods」「Implementation of scTransMIL variants」 | ✅ |
| 07 | 泛癌基准与消融（`07_泛癌基准与消融.md`） | Fig. 2a–c；Supplementary Figs. 1–3 | ✅ |
| 08 | 癌种预测与外部 BRCA 验证（`08_癌种预测与外部BRCA验证.md`） | Fig. 3；Supplementary Figs. 4–5 | ✅ |
| 09 | 单细胞可解释性与基因注意力（`09_单细胞可解释性与基因注意力.md`） | Fig. 2d–f、Fig. 4；Supplementary Figs. 6–7 | ✅ |
| 10 | AML 微调、亚型与生存（`10_AML微调_亚型与生存.md`） | Fig. 5；Supplementary Figs. 8–12 | ✅ |
| 11 | 潜在空间、scIB 与通路富集（`11_潜在空间_scIB与通路富集.md`） | Fig. 3g–j、Fig. 4h | ✅ |
| 12 | 计算性能与可复现性（`12_计算性能与可复现性.md`） | Fig. 2b；Methods「Training」「Statistics」 | ✅ |
| 13 | 复现质控清单、排错与批判性评述（`13_复现质控清单、排错与批判性评述.md`，NN-1） | 全文 | ✅ |
| 14 | 背景知识精讲与方法学原理速查（`14_背景知识精讲与方法学原理速查.md`，NN） | 全文 | ✅ |

---

## 2. 目录与命名（本项目的具体化）

- 项目 slug：`scTransMIL_patient_single_cell_MIL`；目录 `20261018_scTransMIL_patient_single_cell_MIL/`。
- 数据源目录：

| 子目录 | 内容 | 来源 |
|---|---|---|
| `raw/pancancer/` | 104 个 h5ad（泛癌图谱，4 M cells） | https://cellatlas.kaist.ac.kr/ecosystem/（Kang et al. 2024 Nat Commun） |
| `raw/thyroid_GSE184362/` | 甲状腺癌 scRNA（GSE184362） | GEO |
| `raw/brca_cellxgene/` | 外部 BRCA（40 万细胞 / 37 samples） | CELLxGENE `34f5307e-7b4d-4a48-b68f-2ba844c6414b` |
| `raw/aml_GSE116256/` | AML scRNA（35 tumor + 6 normal） | GEO GSE116256 |
| `raw/bulk_aml_dbgap/` | AML bulk RNA-seq（671 例，生存） | dbGaP **phs001657.v2.p1（受控访问）** |
| `raw/pretrained/` | Geneformer `gf-6L-30M-i2048`、token/median 字典、scGPT whole-human 权重 | HuggingFace / 作者 Google Drive |
| `processed/` | bag 索引、tokenized 数据、embedding、splits | 本项目脚本 |
| `results/` | metrics、attention、亚型、生存表 | 本项目脚本 |
| `figures/` | `.pdf` + `.png`（300 dpi） | 本项目脚本 |
| `code/` | `NN_MM_<verb>.py/.sh` | 本项目脚本 |
| `logs/` | pip freeze、运行日志、DECISIONS、review record | 本项目脚本 |

- 命名特例：
  - `slide_id` 的格式沿用作者代码：`{patient}_{tissue}_{bag_index}`（D5）。
  - 标签字符串：`tumor_tissue` / `normal_tissue`（二分类）；癌种用作者 `label_dict` 的 29 个缩写（D16）。
  - 图号前缀 `Fig<N><L>_`；补充图 `FigS<N><L>_`。

---

## 3. 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| `sample` / `patient` / `donor` 混用 | **sample**（一次建库/一个 slide）、**donor**（患者） | MIL 的 bag 是 sample 或 sub-sample；split 以 donor 为单位 |
| `bag` / `sub-sample` 混用 | **bag**（模型输入，≤256 cells）；**sub-sample**（数据构造层） | 原文与代码一致 |
| `tumor score` / `attention` / `A` 混用 | **tumor score**（实例级标量）、**attention weight**（softmax 后权重） | 03 号 |
| `CLAM_SB` / `scTransMIL` | **scTransMIL**（论文名）；`CLAM_SB` 是实现骨架 | 代码 |
| `29 classes` / `30 cancer types` | 正文 **30**；代码 `label_dict` **29** | D16 |
| `n=7` / `k=10` | 正文 **n=7 次重复**；代码默认 `k=10` | D7 |

---

## 4. 产物契约表

| 产物路径 | 生产者(唯一) | 格式/结构（列名·尺度·单位·坐标系） | 消费方 | 生成脚本 |
|---|---|---|---|---|
| `raw/pancancer/*.h5ad` | 外部（01） | AnnData；`obs` 含 `Patient`、`Tissue`、`cnv_status`、`Cancer type`；`var` 为基因（SYMBOL 或 Ensembl，需 D2 统一） | 02 | `code/01_02_fetch_pancancer.sh` |
| `raw/thyroid_GSE184362/*.h5ad` | 外部（01） | 同上的 AnnData；`obs` 含样本/患者/细胞类型 | 09 | `code/01_04_fetch_others.sh` |
| `raw/brca_cellxgene/*.h5ad` | 外部（01） | AnnData；`obs` 含 `sample`、`cell_type`、malignant 注释 | 08 | `code/01_04_fetch_others.sh` |
| `raw/aml_GSE116256/*.h5ad` | 外部（01） | AnnData；`obs` 含 `Sample`、`Patient`、细胞类型、malignant 注释 | 10 | `code/01_04_fetch_others.sh` |
| `raw/bulk_aml_dbgap/*.tsv` | 外部（01，受控） | 表达矩阵 + 生存表（`patient_id`、`time`、`event`） | 10 | `code/01_04_fetch_others.sh` |
| `raw/pretrained/` | 外部（01） | Geneformer `gf-6L-30M-i2048` 权重、`token_dictionary_gc95M.pkl`、`gene_median_dictionary_gc95M.pkl`、scGPT whole-human 权重 | 03、06 | `code/01_01_clone_code.sh` + 手动下载 |
| `processed/gene_id_map.tsv` | 02 | 列：`gene_name`、`ensembl_id`、`token_id`、`mapped` | 02、03 | `code/02_01_make_bags.py` |
| `processed/bag_index_<dataset>.csv` | 02 | 列：`case_id`（donor）、`slide_id`（bag，`{patient}_{tissue}_{idx}`）、`label`、`n_cells`、`dataset`；每个 bag ≤256 cells；由 104 个 part 合并（`02_05_merge_index.py`） | 04、05、07 | `code/02_05_merge_index.py` |
| `processed/bags/<dataset>/<slide_id>/bag.h5ad` | 02 | AnnData：≤256 cells × 全部基因；`obs` 保留 `cnv_status`/`Cancer type`/`Celltype` | 02、03 | `code/02_01_make_bags.py` |
| `processed/tokenized/<dataset>/<slide_id>.pt` | 02 | dict：`input_ids`（int64，N×1024，右 padding）、`attention_mask`（N×1024，0/1） | 03、04 | `code/02_02_tokenize.py` |
| `processed/splits/splits_<task>_<i>.csv` | 05 | **宽表**：`train`/`val`/`test` 三列，每列是该 split 的 `slide_id` 列表（与作者 `save_splits` 格式一致）；按 donor 分层，同一 donor 的所有 bag 只出现在一列 | 04、05 | `code/05_01_make_splits.py` |
| `processed/splits/split_assignments_<task>_<i>.csv` | 05 | 长表：`case_id`、`split`（train/val/test）；用于审计 donor 不重叠 | 05、13 | `code/05_01_make_splits.py` |
| `results/metrics_<task>_<method>_fold<i>.csv` | 05 | 列：`slide_id`、`Y`、`Y_hat`、`Y_prob_*`；二分类/多分类分别存 | 05、07、08 | `code/05_02_eval.py` |
| `results/metrics_summary_<task>.csv` | 05 | 列：`method`、`n_splits`、`acc_mean`、`acc_sd`、`auc_mean`、`auc_sd`、`f1_mean`、`f1_sd` | 07、08、12 | `code/05_03_summarize.py` |
| `processed/embeddings/<dataset>_<split>.npy` | 05 | `float32`，shape = (n_bags, 256)；行顺序与 `results/bag_embedding_index_<dataset>.csv` 的 `slide_id` 一一对应 | 08、11 | `code/05_05_export_embeddings.py` |
| `results/cell_scores_<dataset>.tsv` | 05 | 列：`cell_id`、`slide_id`、`tumor_score`、`cnv_status`、`cell_type` | 09、10 | `code/05_04_infer_cells.py` |
| `results/attention_<dataset>_<cellgroup>.tsv` | 09 | 列：`gene`、`attention_mean`、`rank`（每个细胞组一行/基因） | 09、11、10 | `code/09_03_attention.py` |
| `results/aml_subtypes.csv` | 10 | 列：`cell_id`、`tumor_score`、`subtype`（HSC/GMP/Prog/ProMono/Mono/cDC-like） | 10 | `code/10_03_subtypes.py` |
| `results/survival_aml.csv` | 10 | 列：`patient_id`、`subtype`、`time`、`event`、`risk_group` | 10 | `code/10_05_survival.py` |
| `figures/Fig<N><L>_<slug>.pdf/.png` | 对应模块 | 矢量 PDF + 300 dpi PNG；图注写 n、指标、重复次数 | 00、13 | 各模块脚本 |
| `logs/00_checksums.md5` | 01 | `md5  filename`（raw/ 全量，受控数据只记文件名+大小） | 13 | `code/00_checksum.sh` |

> 纪律（L0 §4.2）：生产者唯一；格式即契约；改格式/名 = 改本表 + 通知消费方 + 记 `DECISIONS.md`。

---

## 5. 术语表（唯一口径）

| 术语 | 本项目定义 |
|---|---|
| **bag** | 模型的一次输入：一个 sample 或一个 ≤256 细胞的 sub-sample；有 bag 级标签 |
| **instance** | bag 内的一个细胞 |
| **MIL 假设** | bag 为阳性 ⇔ 至少一个 instance 阳性；bag 为阴性 ⇔ 所有 instance 阴性 |
| **sub-sample** | 数据构造层把原始 sample 按 256 细胞切块得到；最后一块不足 256 时由代码**回绕复制**补齐（D5） |
| **tumor score（`a_i`）** | 门控分支输出的实例级标量（softmax **前**，正文 Eq.2）；用于 adaptive loss 的排序与伪标签 |
| **attention weight** | 对 `a_i` 在 bag 内做 softmax 后的权重；用于 bag 聚合 `M = A @ h` |
| **cell score** | 用于单细胞 AUC/UMAP 的分数；作者 `eval_utils.py` 保存的是 **softmax 后**的 attention weight（D22）。与 `a_i` 不是同一个量 |
| **attention weight** | `softmax(A)` 后的实例权重；bag 聚合实际使用的是它（D3） |
| **rank value encoding** | 每个细胞去掉零表达基因、按表达量降序排列、映射到 token 的序列；不含绝对表达值 |
| **token/词表** | Geneformer `token_dictionary_gc95M.pkl` 的 Ensembl→ID 映射；`input_ids` 截断到 1024（D4） |
| **CNV status** | 泛癌图谱提供的 corrected CNV 推断标签（`tumor`/`normal`），作为 cell-level 参考真值 |
| **Sample setting** | bag 标签由 bag 内细胞的 CNV 状态决定（默认、低噪声） |
| **Donor setting** | bag 标签由 donor 级诊断决定（高噪声、模拟临床） |
| **split** | 按 **donor** 分层的一次 80/10/10 划分；正文 n=7 次重复（D7） |
| **29 / 30 classes** | 正文 30 cancer types；代码 `label_dict` 29 类（D16） |
| **HSC-like / GMP-like / Mono-like** | AML 亚型：按 tumor score>0.2 筛出恶性细胞后，用 kNN 对齐正常分化状态命名 |
| **scIB** | 单细胞整合基准；本项目把 sample embedding 当 cell、cancer type 当 cell type、dataset 当 batch |

---

## 6. 环境与版本钉死

| 组件 | 版本 | 出处 | 偏离说明 |
|---|---|---|---|
| Python | 3.10 | 作者 README | — |
| PyTorch | ≥2.0（CUDA 11.8/12.x） | 我方选择 | 原文未给版本 |
| transformers | ≥4.30 | Geneformer 依赖 | — |
| Geneformer | HuggingFace `ctheodoris/Geneformer` | 作者 README | 需 `pip install .` |
| 预训练权重 | `gf-6L-30M-i2048` | 作者代码 `./models/gf-6L-30M-i2048/` | 需从作者 Google Drive/Zenodo 获取 |
| token/median 字典 | `*_gc95M.pkl` | 作者代码 | 与 6L-30M 模型的搭配需核对（D10） |
| scGPT（变体） | whole-human 预训练权重 | 原文 Methods | 仅 scTransMIL(FM) 用 |
| scvi-tools | ≥0.20 | 基线 scVI+MIL | — |
| smooth-topk | GitHub `oval-group/smooth-topk` | 作者 README | 自适应聚类损失 |
| scanpy/anndata | ≥1.9 / ≥0.10 | 预处理 | — |
| scIB | ≥1.1 | 整合评估 | — |
| GSEApy | ≥1.1 | 通路富集 | — |
| lifelines / scikit-survival | 最新 | 生存分析 | 原文未指定 |
| R | 4.3+（仅 DEGAS 基线） | 作者用原 R 包 | DEGAS 是 R 包 |
| GPU | NVIDIA TESLA V100（32 GB） | 原文 Fig. 2b | 显存不足时降 batch/序列长度并声明 |
| CPU/RAM | 10 cores / 40 GB | 作者 README | — |

> 每脚本结束写 `logs/<NN>_pip_freeze.txt` 与 `logs/<NN>_gpu.txt`（`nvidia-smi`），并记录完整执行命令（L0 §5.6）。

---

## 7. 统计与随机性

- **随机种子常量**：`20261018`
  - Python：`random.seed(20261018)`、`np.random.seed(20261018)`、`torch.manual_seed(20261018)`、`torch.cuda.manual_seed_all(20261018)`。
  - 作者代码 `seed_torch(seed=1)`；本项目在**每次 split** 用 `20261018 + fold`，并记录（D18）。
- 必须固定种子的环节：donor 划分、bag 采样、weighted sampler、dropout、模型初始化、UMAP、kNN、生存分析的分层。
- 项目特有统计约定：
  1. **主指标**：样本级 Accuracy / AUC / macro-F1，**mean ± s.d. across n=7 splits**（与原文一致）。
  2. **单细胞指标**：以 CNV status 为参照的 ROC-AUC；报告 Acc/F1 时注明阈值（0.5）。
  3. **多分类**：macro-F1；混淆矩阵与 Sankey 图并报。
  4. **不确定性**：所有均值必须带 s.d. 与 n；不得只报最好一次 split。
  5. **生存**：KM + log-rank；报告 n、事件数、中位生存；Cox 需 PH 检验（若使用）。
  6. **富集**：Fisher 单侧 + BH；背景基因集必须与词表一致（D15）。

---

## 8. 决策登记表（D 表索引）

> **权威全文见 `13_复现质控清单、排错与批判性评述.md` §八**；本表只列编号与一句话裁决。

| # | 决策项 | 涉及模块 | 裁决(主路径) | 备选/敏感性 | 依据 | ★ |
|---|---|---|---|---|---|---|
| D1 | 泛癌数据版本与获取 | 01、02、07 | 用 `cellatlas.kaist.ac.kr/ecosystem` 的 Kang et al. 2024 图谱；按作者 104 个 h5ad 清单核对 | 用 CELLxGENE 镜像 | 原文 Data availability；作者 README | ★ |
| D2 | 基因 ID 与词表 | 02、03 | 统一 Ensembl；用作者 `gene_id.csv` 把 SYMBOL 映射到 Ensembl，再映射到 `token_dictionary_gc95M.pkl`；未映射基因丢弃并记录 | 用新版 Ensembl | 作者 `1-process_v3.py` | ★ |
| D3 | 正文 vs 代码的架构差异 | 03、04 | **以代码为准复现**：BERT `pooler_output`（非正文的 mean pooling）；bag 聚合 `softmax(A) @ h`（非正文的 `Σ a_i h_gate`）；门控 `L=256, D=64`；两阶段为两个独立任务；`a_i`=softmax 前、`cell_score`=softmax 后（D22）；instance 分类器输入原始 `h`（非 `h_gate`） | 按正文公式另写一版做对照 | 作者 `models/model_clam.py`；正文 Eq.1–6 | ★ |
| D4 | 序列长度 | 02、03 | `max_len=1024`（代码截断）；正文 "full transcriptome" 是相对 PCA 变体而言 | 2048（Geneformer 30M 的 model_input_size） | 作者 `1-process_v3.py`；Geneformer 文档 | ★ |
| D5 | sub-sample 大小与尾部补齐 | 02 | 256 cells/bag；最后一块不足 256 时**回绕复制**开头细胞补齐；<17 cells 的 bag 在 step2 过滤 | 丢弃最后一块 | 作者 `split_ids()`、step2 | ★ |
| D6 | CNV 标签来源 | 02、07、09 | 直接用图谱提供的 corrected CNV status；不重跑 inferCNV | 用 inferCNVpy 重算做敏感性 | 原文 Methods；图谱 | ★ |
| D7 | 数据划分 | 05 | 按 donor 分层 80/10/10；**n=7 次独立随机划分**（正文） | 代码默认 10 folds | 原文 Methods；代码 `k=10` | ★ |
| D8 | 类别不平衡 | 04 | weighted sampler：采样概率 ∝ 1/类别频率 | 不做加权 | 原文 Methods；`--weighted_sample` | |
| D9 | 基线版本 | 06 | KIDA/scIDST 用官方 Python；DEGAS 用官方 R 包；scVI+MIL 用 scvi-tools + MultiMIL attention | 用论文报告值替代 | 原文 Methods；Supplementary Data 5 | ★ |
| D10 | 预训练权重与字典 | 01、03 | Geneformer `gf-6L-30M-i2048` + `gc95M` 字典（按作者代码）；记录 md5。**该 checkpoint 在当前 HF 仓库未直接列出**（现有 V1-10M / V2-104M / V2-316M），需从 V1 归档/Zenodo/作者 Drive 获取 | 用 `Geneformer-V1-10M`（同为 6 层/hidden 256）替代 | 作者 README/代码；HF 文件列表；`原文未报告`搭配理由 | ★ |
| D11 | 硬件与计时口径 | 12 | V100 32GB、10 cores、40GB；计时含预处理/特征提取；报相对值 | 其他 GPU | 原文 Fig. 2b；README | |
| D12 | bulk AML 受控数据 | 10 | dbGaP `phs001657.v2.p1` 需申请；不可得则生存分析标为"未复现" | 用 GEO 其他 AML 队列 | 原文 Data availability | ★ |
| D13 | AML 微调策略 | 10 | **冻结 embeddings + 前 4 层，训练最后 2 个 Transformer 块 + pooler + aggregator**，30 epochs（与 `finetune_utils.py` 一致）；正文 Methods 的"freeze instance encoder"与 Results 的"last two blocks"冲突，按代码落地 | 只训 aggregator；全量微调 | 正文 pp.10/14；`utils/finetune_utils.py:174-187`（`def_freeze_layers=4`） | ★ |
| D14 | AML 亚型与生存 | 10 | tumor score>0.2 → kNN 对齐正常分化；每亚型取 top 50 患者（top 30 基因均表达）→ 6×50=300；KM | 不同阈值/患者数 | 原文 Methods/Results | ★ |
| D15 | 富集工具与背景 | 09、11 | GSEApy `enrichr`，Fisher 单侧，BH；背景用词表内基因 | 超几何 + 自定义背景 | 原文 Methods | |
| D16 | 29 vs 30 类 | 05、08 | 代码 `label_dict` 为 **29 个癌种标签**；`2-bag_eval.py` 把 Stage-1 判 normal 的 bag 置为 `Y_hat=29.0`，形成第 30 个**预测桶（normal/rejection）**，不一定是"缺失的癌种"。报告时明确该桶含义；若需 30 个癌种，先用 Supplementary Data 1 核对 | 强行 30 类 | 代码 `main.py:171`、`2-bag_eval.py:24-30`；正文 30 | ★ |
| D17 | 细胞 QC | 02 | 沿用图谱已 QC 的 h5ad；不额外过滤（除零表达基因） | 统一重做 QC | 原文用公开 processed 矩阵 | |
| D18 | 随机种子 | 全部 | `20261018 + fold`；每次 split 记录 | 作者 `seed=1` | 本项目规范 | |
| D19 | 训练超参 | 04 | 以正文为准：lr 2e-4、weight decay 1e-5、30–200 epochs、patience 10；覆盖 `main.py` 默认（5e-5/30） | 用代码默认 | 正文 Methods；代码 `main.py` | ★ |
| D20 | AUC 口径 | 05 | 用 `Y_prob`（软概率）算 AUC | 沿用 `2-bag_eval.py` 的 `Y_hat` | 原文 Methods；代码 `2-bag_eval.py` | ★ |
| D21 | 基因 attention 归约 | 09 | 论文只给"最后一层 + 所有 head 平均 + top-100"；**行/列归约方式论文与代码均未给出**（`model_clam.py` 硬编码 `output_attentions=False`，仓库无提取脚本）。本项目**假设**取 CLS→token 的注意力，并做"列求和"敏感性对照 | 列求和/其他行 | 正文 Methods「Application」；代码无实现 | ★ |
| D22 | cell score 用哪个量 | 05、09、10 | 单细胞 AUC/UMAP 用作者 `eval_utils.py` 保存的 **softmax 后 attention weight**（这是 0.71 的来源）；adaptive loss 用 softmax 前 `a_i`；两者都输出、分别命名 | 用 softmax 前 `A_raw` 重算 AUC | `utils/eval_utils.py:87-91`；正文 Eq.2 | ★ |
| D23 | 甲状腺数据是否在训练集 | 01、09 | 作者 README 的泛癌清单包含 `tc_ptc_GSE184362.h5ad` → **Patient 5 可能已在训练集**。主路径训练时排除 GSE184362，再做"不排除"的对照；若无法排除，09 的 0.64 标为 in-corpus case study | 保持原训练集并声明 | 作者 README 文件清单；正文称 external validation | ★ |

---

## 9. QC 门槛

- **数值层**：原文正文给出的数字（4 M/1,563/999/30/104/30k、0.71/0.67/0.58/0.56、245.16 s/720 s、0.64、35+6/30k/21 donors、671 例）必须能对上；图上的 Acc/AUC/F1 以 Source Data 为准（正文未逐项给数）。
- **判据层**：✅ 通过 → 推进；⚠️ 小幅偏差 → 排查后记录放行；❌ 显著偏差 → **第二套独立实现**（如按正文公式重写聚合）复核，仍不符则入 gap 表。
- **结题层**：gap 表逐模块填「文献值 / 我方复现值 / 是否一致 / 差异原因 / 处理方式」。
- 全局 Go/No-Go（摘要）：
  - G1 数据完整性：泛癌 4 M cells / 1,563 samples / 999 donors / 30 cancer types / 104 datasets（±5%）。
  - G2 样本级性能：scTransMIL 在 Acc/AUC/F1 上**优于或持平**所有基线（相对判据，不依赖具体绝对值）。
  - G3 单细胞：AUC ≥0.68（原文 0.71，允许 −0.03）；高于 scVI+MIL（0.67）。
  - G4 甲状腺：Patient 5 paratumor 的预测分数在 0.55–0.75（原文 0.64），且高 score 富集在上皮细胞。
  - G5 AML：6 个亚型可复现；HSC-like 的生存显著差于 GMP-like（log-rank p<0.05）。
  - G6 效率：单 bag 推理时间与 DEGAS 的比值 <0.5（原文 245/720 ≈0.34）。

---

## 10. 偏离声明（对 L0 的非红线偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| §8.1 头部四要素 | 00/13/14 号不强制"模块对应原文"四要素 | L0 §1.3 将 00/NN-1/NN 定义为非模块文档 | 仅这三篇 |
| §6.4 下载校验 | 泛癌 104 个 h5ad 体量巨大，只对实际使用的文件做 md5；dbGaP 受控数据只记文件名/大小 | 存储与访问限制 | `logs/00_checksums.md5` 头部注明 |
| §3.2 图像双份 | 中间调试图可只出 PNG；最终图必须 PDF+PNG | 训练调试图数量多 | `figures/` 最终图不受影响 |

> L0 §12 诚实性红线不可偏离。

---

## 11. 修订记录

| 版本 | 日期 | 修改内容 | 修改人 |
|---|---|---|---|
| v1.0 | 2026-10-18 | 首版：15 篇文档（00–14）、产物契约表、术语表、D1–D23、环境版本表 | 复现项目组 |
