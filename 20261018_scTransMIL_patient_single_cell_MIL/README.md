# scTransMIL 患者级-单细胞多实例学习 —— 文献复现项目

**文献**：Tang Z, Wang F, Yang F, Song J, Li Y, Zhou J, Song Y, Chen S, Zhu J, You L, Chen CY-C, Yao J. *scTransMIL bridges patient-level disease states and single-cell transcriptomics for cancer screening and heterogeneity inference.* **Nature Communications (2026) 17:10800.** DOI: [10.1038/s41467-026-77538-5](https://doi.org/10.1038/s41467-026-77538-5)

**一句话**：以作者公开代码 `TencentAILabHealthcare/scTransMIL` 为基准，重训并核验"rank 编码 + 6 层 Transformer + 门控注意力 + 多实例学习"流水线，在泛癌、甲状腺、外部 BRCA、AML 四个数据域上复现样本级诊断、单细胞 tumor score、癌种分类、注意力基因与 AML 亚型/生存，并逐条登记正文与代码不一致处（D1–D23）。

---

## 一、文件索引（本目录）

| 文件 | 类型 | 一句话描述 |
|---|---|---|
| `README.md` | 索引 | 本文件：项目一句话、文件清单、使用顺序 |
| `PROJECT_SPEC.md` | **L1 规范** | 模块表、产物契约表、术语表、D1–D23 索引、环境版本、偏离声明 |
| `00_目录与复现总览.md` | 总览 | 文档目录、科学问题、证据链 C1–C10、全流程图、图表映射、数据资产、算力工期、全局门槛 |
| `01_环境搭建与数据获取.md` | 模块 | Python 3.10/PyTorch/Geneformer 环境；泛癌/甲状腺/BRCA/AML/dbGaP/预训练权重；校验和 |
| `02_数据预处理与MIL样本构造.md` | 模块 | gene_id 映射、rank 编码、256-cell bag（尾部回绕）、tokenize 1024、bag 索引与标签 |
| `03_模型架构_rank编码_Transformer_门控与聚合.md` | 模块 | 6 层 Transformer、门控 `tanh⊙sigmoid`、tumor score、softmax 聚合；**正文 vs 代码差异** |
| `04_损失函数与训练策略.md` | 模块 | bag BCE/CE、adaptive clustering loss（B=8、τ=1、α=1、β=0.7）、SmoothTop1SVM、训练超参 |
| `05_推理_评估与数据划分.md` | 模块 | sub-sample→sample 聚合、多分类投票、Acc/AUC/macro-F1、donor 80/10/10、n=7、单细胞 ROC |
| `06_基线方法与变体复现.md` | 模块 | KIDA、DEGAS、scIDST、scVI+MIL、scTransMIL(PCA/FM) |
| `07_泛癌基准与消融.md` | 模块 | Fig. 2a–c：训练量、推理时间、三个消融、gating 指标、Donor setting |
| `08_癌种预测与外部BRCA验证.md` | 模块 | 30 类原发/转移 Sankey、BRCA 37 samples、误分类、bag UMAP |
| `09_单细胞可解释性与基因注意力.md` | 模块 | 单细胞 AUC 0.71、甲状腺 Patient 5（0.64）、ALL B 细胞、BLK、GSEApy |
| `10_AML微调_亚型与生存.md` | 模块 | GSE116256 微调、6 亚型、纵向治疗、bulk 671 生存 |
| `11_潜在空间_scIB与通路富集.md` | 模块 | bag UMAP、scIB 八项、富集 top-20 |
| `12_计算性能与可复现性.md` | 模块 | 245.16 s vs 720 s、显存、7 次重复、种子 |
| `13_复现质控清单、排错与批判性评述.md` | **L2 + D 表** | N1–N25 门槛、G1–G6、gap 表、排错、**D1–D23 权威全文**、证据等级与批判 |
| `14_背景知识精讲与方法学原理速查.md` | 背景 | MIL 谱系、Transformer、rank 编码、SmoothTop1SVM、scIB、AML 生物学、FAQ |
| `2026_scTransMIL_..._inference.pdf` | 原文 | 文献 PDF（17 页） |
| `logs/DECISIONS.md` | 日志 | 运行时决策流水 |
| `logs/13_review_record.md` | 审核记录 | S4 独立审核的发现→修正对照，及待 S6 裁决的 8 项 |

## 二、目录结构

```text
20261018_scTransMIL_patient_single_cell_MIL/
├── README.md                  # 本索引
├── PROJECT_SPEC.md            # L1 项目规范
├── 00 … 14                    # 15 篇文档
├── 2026_scTransMIL_...pdf     # 原文
├── raw/
│   ├── pancancer/             # 104 个 h5ad
│   ├── thyroid_GSE184362/
│   ├── brca_cellxgene/
│   ├── aml_GSE116256/
│   ├── bulk_aml_dbgap/        # dbGaP 受控
│   └── pretrained/            # Geneformer/scGPT/作者权重
├── processed/                 # bags/tokenized/splits/embeddings
├── results/                   # metrics/attention/subtypes/survival
├── figures/                   # .pdf + .png
├── code/                      # NN_MM_*.py/.sh（S6 实现）
└── logs/                      # pip freeze/nvidia-smi/DECISIONS/review
```

## 三、使用顺序

1. **先读** `13` 号 §八（D1–D23）与 §九（批判）→ 知道哪些是正文没写、由我们裁决的，尤其是 **D3（正文公式 vs 代码实现）**。
2. **再读** `14` 号（MIL/Transformer/rank 编码/损失速查）。
3. **执行**：`01 → 02 → 03 → 04 → 05` 是主链；`06` 基线；`07/08` 评估；`09/10` 应用；`11/12` 表征/性能。
4. **结题**：按 `13` 号 §五 填 QC、§七 填 gap、§十 打钩。

## 四、状态

| 阶段 | 状态 |
|---|---|
| S1 立项 / S2 规范 / S3 文档 | ✅ 完成（15 篇） |
| S4 独立审核 | ✅ 已完成（2 CRITICAL / 17 HIGH 已修正；记录见 `logs/13_review_record.md`） |
| S5 自检 | ✅ `tools/repro_check.sh` 全项通过 |
| S6 执行（下载数据、训练、填 QC/gap） | ⏳ 未开始 |

> **诚实性提醒**：正文未逐项给出 Fig. 2/3 的 Acc/AUC/F1 数值（只有 Source Data）；本项目不编造这些数。凡标"原文未报告"处，均为复现者裁决并登记在 D 表。不要为了凑数改参数（L0 §12）。
