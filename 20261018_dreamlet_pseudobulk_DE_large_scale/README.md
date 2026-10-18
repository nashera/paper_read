# Dreamlet 大规模单细胞差异表达 —— 文献复现项目

**文献**：Hoffman GE, Lee D, Bendl J, … Roussos P. *Efficient differential expression analysis of large-scale single-cell transcriptomics data using Dreamlet.* **Nature Communications (2026) 17:9597.** DOI: [10.1038/s41467-026-75680-8](https://doi.org/10.1038/s41467-026-75680-8)

**一句话**：用作者公开的 `dreamlet` R 包与 `dreamlet_analysis` 代码，重跑"pseudobulk + 两阶段精度权重 + 线性混合模型 + 经验贝叶斯 moderated t"流水线，在 COVID-19、前列腺癌骨转移、阿尔茨海默病三个公开队列上逐图核验原文结论，并补齐原文未写清的参数与公式（D1–D20）。

---

## 一、文件索引（本目录）

| 文件 | 类型 | 一句话描述 |
|---|---|---|
| `README.md` | 索引 | 本文件：项目一句话、文件清单、使用顺序 |
| `PROJECT_SPEC.md` | **L1 规范** | 模块编号表、产物契约表、术语表、D 表索引（D1–D20）、环境版本、偏离声明 |
| `00_目录与复现总览.md` | 总览 | 文档目录、科学问题、证据链 C1–C10、全流程图、图表→文档映射、数据资产、算力工期、全局 QC 门槛 |
| `01_环境搭建与数据获取.md` | 模块 | R 4.6/Bioc 3.23/dreamlet 1.10.0 安装；COMBAT、GSE143791、PsychAD（syn51188606）的获取与 h5ad 体检；md5 校验 |
| `02_数据接入与pseudobulk聚合.md` | 模块 | H5AD→SCE（zellkonverter/DelayedArray）、`aggregateToPseudoBulk()` 的分块/beachmat 原理、内存量级推导、细胞数矩阵 |
| `03_两阶段精度权重与均值方差趋势.md` | 模块 | Poisson 初始化权重的 δ 法推导（含方差膨胀因子）、TMM + 自适应先验计数、voom 趋势与 `1/ŝ⁴`、两级 vs 单级对照 |
| `04_精度加权LMM与经验贝叶斯moderated_t.md` | 模块 | LMM/GLS/REML、为何不用固定效应/GLMM、残差自由度 `n−tr(H)`、Smyth EB（d0、s0²、moderated t）、模拟单元测试 |
| `05_差异表达分析流程与结果提取.md` | 模块 | `processAssays()`→`dreamlet()`→`topTable()`；两个公式（权重 vs 检验）；contrasts；细胞类型内 BH vs study-wide BH；π1；火山/森林/热图 |
| `06_方差分解与批次效应量化.md` | 模块 | `fitVarPart()` 原理与公式；COVID/Kfoury/AD 的 VP 公式；代表基因；GC 含量与 pool 方差分数的 Spearman 相关 |
| `07_基因集分析camera与zenith.md` | 模块 | camera 竞争性检验（VIF 校正）、zenith 对 LMM 的适配、GO/MSigDB、`plotZenithResults()`、手算复现 |
| `08_计算与统计性能基准复现.md` | 模块 | Fig. 2A–E：5 种 pseudobulk 方法、6 种 DE 方法、muscat 模拟（AUPR/F1）、真实数据置换 FPR |
| `09_案例I_COVID19重症单核细胞.md` | 模块 | COMBAT 674K/110：样本筛选（最重症 + Oxford）、cMono VP、1824 DEG、CLU/NFKBIA、TNFα/NF-κB |
| `10_案例II_前列腺癌骨转移.md` | 模块 | Kfoury GSE143791：重复测量 LMM、三个 contrast、157 DEG、MHC-II/chaperone |
| `11_案例III_AD队列数据生成与上游处理.md` | 模块 | 湿实验 SOP（蔗糖梯度/FANS/hashing/10x）+ 干实验（STARsolo/cellSNP+vireo/QTLtools-mbv/QC/Harmony/Leiden/注释/sex check） |
| `12_案例IV_AD差异表达与下游分析.md` | 模块 | PsychAD 586/299/1.4M/22 类：VP、1037 DEG、PTPRG 1.59、PDE10A、DUSP10、p38MAPK、核数-效力 |
| `13_方法对照与跨队列一致性.md` | 模块 | dreamlet vs limma vs DESeq2；与 Mathys 2023 跨队列；10X pool 随机效应 vs Harmony 前 5 成分 |
| `14_复现质控清单、排错与批判性评述.md` | **L2 + D 表** | 26 项数值门槛、G1–G7 判据、gap 表、跨模块排错、**D1–D20 权威全文**、证据等级与批判性评述 |
| `15_背景知识精讲与方法学原理速查.md` | 背景 | 方法谱系、LMM/voom/EB/VP/camera/BH 的数学速查、AD/COVID/骨转移生物学、FAQ、文献导读 |
| `2026_Efficient_..._Dreamlet.pdf` | 原文 | 文献 PDF（12 页） |
| `logs/DECISIONS.md` | 日志 | 运行时决策流水（L0 §10.2） |
| `logs/14_review_record.md` | 审核记录 | S4 独立审核的发现→修正对照，及待 S6 实测裁决的 6 项 |

## 二、目录结构

```text
20261018_dreamlet_pseudobulk_DE_large_scale/
├── README.md                  # 本索引
├── PROJECT_SPEC.md            # L1 项目规范
├── 00 … 15                    # 16 篇文档（总览 / 模块 / 质控 / 背景）
├── 2026_..._Dreamlet.pdf      # 原文
├── raw/                       # 只读原始数据
│   ├── combat_covid19/        # COMBAT h5ad + 临床变量
│   ├── kfoury_prostate/       # GSE143791 h5ad/10x
│   ├── psychad_ad/            # PsychAD h5ad（syn51188606）
│   └── benchmark/             # 基准 H5AD、MSigDB/GO、GC 参考
├── processed/                 # 中间对象（pb_*.rds、proc_*.rds、fit_*.rds、vp_*.rds）
├── results/                   # 结果表（topTable_*.tsv、vp_*.tsv、zenith_*.tsv、benchmark_*.csv）
├── figures/                   # 图（.pdf + .png，300 dpi）
├── code/                      # 脚本（NN_MM_*.R/.py/.sh；S6 实现）
└── logs/                      # sessionInfo、运行日志、DECISIONS.md、md5
```

## 三、使用顺序

1. **先读** `14` 号 §八（D1–D20）与 §九（批判性评述）→ 知道哪些参数是原文没写、由我们裁决的。
2. **再读** `15` 号（方法学速查）→ 补齐 LMM/voom/EB/VP/camera 的数学背景。
3. **执行**：`01 → 02 → 03 → 05 → 06 → 07` 是通用管线；`09/10/12` 是三个案例；`08` 是性能基准；`11` 是 AD 上游（可选重跑）；`13` 是对照。
4. **结题**：按 `14` 号 §五 填 QC 表、§七 填 gap 表、§十 打钩。

## 四、状态

| 阶段 | 状态 |
|---|---|
| S1 立项 / S2 规范 / S3 文档 | ✅ 完成（16 篇） |
| S4 独立审核 | ✅ 已完成（code-reviewer + verifier 双通道；8 CRITICAL/18 HIGH 已修正，记录见 `logs/14_review_record.md`） |
| S5 自检 | ✅ `tools/repro_check.sh` 通过 |
| S6 执行（下载数据、跑代码、填 QC/gap 表） | ⏳ 未开始 |

> **诚实性提醒**：本项目文档中标注"原文报告"的数值均可在 PDF 正文定位（少数参数来自作者代码，已注明）；凡标"原文未报告"处，均为复现者裁决并登记在 D 表。不要为了凑数改参数（L0 §12）。
