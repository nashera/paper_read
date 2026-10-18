# CRC 髓系靶向治疗单细胞复现 —— 文献复现项目

**文献**：Zhang L, Li Z, Skrzypczynska KM, … Egen JG, Zhang Z, Yu X. *Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer.* **Cell. 2020;181(2):442–459.e29.** DOI: [10.1016/j.cell.2020.03.048](https://doi.org/10.1016/j.cell.2020.03.048)

**一句话**：用原文公开的人类 CRC 单细胞数据（GEO GSE146771）、小鼠 anti-CSF1R / anti-CD40 处理数据（ENA PRJEB34105、ArrayExpress E-MTAB-8832）与 TCGA/GTEx bulk 数据，逐图重建"髓系图谱 → C1QC⁺/SPP1⁺ TAM 二分 → 细胞互作网络 → 髓系靶向治疗机制"这条证据链，并补齐原文 Methods 未写清的参数与算法（D1–D22）。

---

## 一、文件索引（本目录）

| 文件 | 类型 | 一句话描述 |
|---|---|---|
| `README.md` | 索引 | 本文件：项目一句话、文件清单、目录结构、使用顺序 |
| `PROJECT_SPEC.md` | **L1 规范** | 模块编号表、产物契约表、术语表、D 表索引（D1–D22）、环境版本、偏离声明 |
| `00_目录与复现总览.md` | 总览 | 文档目录、科学问题、证据链 C1–C8、全流程图、图表→文档映射、数据资产、算力工期、全局门槛 G1–G13 |
| `01_环境搭建与数据获取.md` | 模块 | R 4.5.1/Bioc 3.22/Python 3.11 环境；GSE146771、PRJEB34105、E-MTAB-8832、TCGA/GTEx 的获取、校验与目录初始化 |
| `02_人类CRC数据接入与统一预处理.md` | 模块 | Smart-seq2 kallisto+tximport、10x Cell Ranger；QC 阈值；熵指标；pagoda2 batch factor；Scrublet + CD14⁺IL32⁺ 双细胞簇 |
| `03_无监督聚类注释与恶性细胞CNV.md` | 模块 | Seurat 图聚类、38/36 免疫簇、13 髓系；跨平台 logistic regression；Ro/e；CNV 推断（moving average + baseline） |
| `04_髓系亚群与TAM二分法.md` | 模块 | hM01–hM13 marker；C1QC⁺ vs SPP1⁺ DEG/GSVA/M1-M2/增殖；多色 IHC |
| `05_单核巨噬细胞发育轨迹.md` | 模块 | RNA velocity（destiny+velocyto）、URD、PAGA；FCN1⁺→TAM；人单核细胞体外分化 |
| `06_细胞互作网络I_scRNA_bulk相关性网络.md` | 模块 | signature → TCGA/GTEx 丰度 → Pearson → self-expressed 置 0 → top13 → z>1.96 网络边 |
| `07_细胞互作网络II_配体受体与置换检验.md` | 模块 | FANTOM5∪CellPhoneDB；E_L×E_R；1,000 次置换；CXCL10–CXCR3、SPP1–ITGAV、SDC2–MMP2、FN1–ITGA5 |
| `08_小鼠肿瘤模型与跨物种髓系整合.md` | 模块 | MC38/Renca 模型与给药；Seurat v3 整合 15 个髓系簇 + T 细胞 mT01–mT12；人–鼠 logistic regression 相似性 |
| `09_antiCSF1R_TAM差异耗竭与耐药.md` | 模块 | mM12/mM14 敏感、mM15 抵抗；Vegfa/Cd274/Arg1↑；血管贴附；TCGA C1QC/SPP1 生存 |
| `10_antiCD40_cDC1激活与CD40_CD40LG轴.md` | 模块 | mM07 cDC1-Ccl22 扩增；CD80/CD86/IL-12↑；activated DC signature 与 OS |
| `11_antiCD40_T细胞动态TCR与STARTRAC.md` | 模块 | Tem↑/Tex↓；TCR αβ；STARTRAC expansion/migration/transition；Bhlhe40⁺ Th1-like 扩增 |
| `12_TCGA_GTEx_bulk验证与生存分析.md` | 模块 | 55/45 分位分组；Cox 校正 stage+MSI；KM；MSI-H；Th1-like–cDC1 相关 |
| `13_实验验证与湿实验SOP.md` | 模块 | 多色 IHC/IF、流式、胞内细胞因子、OTII 培养、单核分化、MSI PCR、外显子测序 |
| `14_复现质控清单、排错与批判性评述.md` | **L2 + D 表** | QC 门槛、gap 表、排错、**D1–D22 权威全文**、证据等级与批判性评述 |
| `15_文献精读、生物学背景与方法学原理速查.md` | 背景 | 逐图精读；CRC 免疫微环境/髓系发育/TAM 二分/cDC1–Th1–CD8；轨迹/网络/LR/STARTRAC/Cox 数学速查 |
| `2020_Zhang_L_Cell_Single-Cell_Analyses_Inform_Mechanisms_of_Myeloid-Targeted_Therapies_in_Colon_Cancer.pdf` | 原文 | 47 页 PDF（正文 + STAR Methods + 补充图 S1–S8） |
| `logs/DECISIONS.md` | 日志 | 运行时决策流水（L0 §10.2） |
| `logs/14_review_record.md` | 审核记录 | S4 独立审核的发现→修复对照、复审结论 |

## 二、目录结构

```text
20261018_CRC_myeloid_scRNA_CSF1R_CD40/
├── README.md                  # 本索引
├── PROJECT_SPEC.md            # L1 项目规范
├── 00 … 15                    # 16 篇文档（总览 / 模块 / 质控 / 背景）
├── 2020_..._Colon_Cancer.pdf  # 原文
├── raw/                       # 只读原始数据
│   ├── human_crc_10x/         # GSE146771 10x processed matrix
│   ├── human_crc_smartseq2/   # GSE146771 Smart-seq2 processed matrix
│   ├── mouse_mc38/            # E-MTAB-8832 / PRJEB34105
│   ├── mouse_renca/           # E-MTAB-8832 / PRJEB34105
│   ├── tcga_coadread/         # UCSC Xena
│   └── gtex_colon/            # UCSC Xena
├── processed/                 # 中间对象（qc/norm/annot/cnv/trajectory）
├── results/                   # 结果表（markers/gsva/network/lr/survival）
├── figures/                   # 图（.pdf + .png，300 dpi）
├── code/                      # 脚本（NN_MM_*.R/.py/.sh；S6 实现）
└── logs/                      # sessionInfo、运行日志、DECISIONS.md、md5
```

## 三、使用顺序

1. **先读** `15` 号（逐图精读 + 生物学背景 + 方法学速查）→ 先建立"这篇文章在讲什么、证据链怎么走"的整体认识。
2. **再读** `14` 号 §十一（D1–D22）与 §十（批判性评述）→ 知道哪些参数是原文没写、由我们裁决的，哪些结论不能升级。
3. **执行**：`01 → 02 → 03 → 04 → 05` 是人类图谱主线；`06 → 07` 是互作网络；`08 → 09 → 10 → 11` 是小鼠机制；`12` 是 bulk/生存；`13` 是湿实验 SOP（本项目不执行）。
4. **结题**：按 `14` 号 §六 填 QC 表、§九 填 gap 表、§十二 打钩；跑 `tools/repro_check.sh` 通过后归档。

## 四、边界声明

- 本项目**不执行湿实验**（D19）：13 号文档给 SOP 与验收标准，gap 表标注"未执行"。
- 人类原始 FASTQ 因隐私不公开；只使用 GEO processed matrix（D1）。
- 所有轨迹/网络/LR/生存结论均为**相关性/预测性**，不得升级为因果（L0 §12）。
- 原文未报告的参数一律显式标注「文献未报告」，并按 D 表裁决；禁止反向调参凑文献数值。
