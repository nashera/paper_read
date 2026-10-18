# 项目决策流水

| 日期 | 决策 | 触发原因 | 依据（D 编号 / 文献 / 实测） | 影响范围 | 决策人 |
|---|---|---|---|---|---|
| 2026-10-18 | 建立项目目录 `20261018_resistance_continuum_cancer_therapy/`，按 L0 §3.1 手工创建 `raw/ processed/ results/ figures/ code/ logs/` | `tools/repro_onboard.sh` 在本机 macOS 上因 `date -d` 不可用而报"不是合法日历日期"，无法执行 | L0 §3.1、附录 B 第 4 步；实测脚本报错 | 全局 | — |
| 2026-10-18 | 原文 PDF 重命名为 `Cellular adaptation to cancer therapy along a resistance continuum (Nature 2024, 631-876).pdf` 并归档到项目根 | L0 §1.2 原文归档；原文件名含中文且与文献内容不匹配 | L0 §1.2、§3.2 第 6 条（PDF 豁免空格） | 全局 | — |
| 2026-10-18 | 数据主路径使用 GEO 处理后 counts/metadata，而非从 SRA 重跑 raw fastq | GSE247691 8 个子系列已公开处理矩阵；raw 路线算力/存储成本高 | D2；GEO 文件盘点（01 号文档 §2.2） | 03–08 | — |
| 2026-10-18 | 记录 ATAC 数据缺口：GSE247685 仅公开 6 个系列级 bigWig + consensus peaks，28 个样本的 per-sample 计数矩阵未公开 | SOFT 中各样本 `supplementary_file_1 = NONE`；系列级 suppl 仅 6 个 bw | D7；实测下载清单 | 06 | — |
| 2026-10-18 | 固定 Seurat 4.4.0（R 4.4.2），而非原文的 4.0.5/2.3.2 | 原文版本组合异常（2.3.2 与 4.0.5 混用），且 4.0.5 在新 R 上兼容性差 | D1；Methods「Single-cell data analysis」 | 03/04/05 | — |
| 2026-10-18 | 状态数主口径 k=5，k=4/6 作敏感性并报告 ARI | 原文明确 k=5，但未给切树方式 | D4；Methods「Defining cell states」 | 04 | — |
| 2026-10-18 | 模块识别按 Neftel 流程重建，目标 6 个模块；候选簇枚举 k=2..20 | 原文只说"recovered all clusters"，未给枚举范围 | D5/D16；Methods「Gene module identification」；作者仓库 `modules/gene_modules.R` | 04 | — |
| 2026-10-18 | 体外群体统一称 adapted line；仅 PDX 三轮后称 resistant | 原文自己也指出最适应系仍对适应剂量敏感 | D11；Results/ Discussion | 全局 | — |
| 2026-10-18 | 作者代码仅作参数对照；缺失脚本按 Methods 重建，不伪造 `.rda` 中间对象 | 仓库脚本 `load()` 本地对象但未发布 | D14；GitHub `yanailab/drug-resistance` 实测 | 01/全部 | — |
| 2026-10-18 | 湿实验模块（02/09）只交付 SOP 与统计模板，不在本项目执行 | 项目定位为文献复现方案与计算复现；无细胞房/动物房条件 | D13；任务范围 | 02/09 | — |
| 2026-10-18 | 代谢组先核验 ST003102 可下载性；不可得则用补充表做定性核对 | Metabolomics Workbench 数据获取方式未确认 | D8；Data availability | 08 | — |
| 2026-10-18 | 全项目随机种子固定为 `20261018` | L0 §6.1 要求统一常量 | D15；L0 §6.1 | 全局 | — |
| 2026-10-18 | 更正 ATAC 数据盘点：GSE247685 consensus 文件实为 chr/start/end + 20 样本计数；Fig. 3c–e 用公开计数直接跑 DESeq2，Fig. 3f 的 8 个 T320/T320W 加药样本走 SRA | 2026-10-18 实测文件表头；此前"per-sample 计数未公开"的判断不准确 | D7；GSE247685 文件实测 | 06 | — |
| 2026-10-18 | 更正单细胞 metadata 口径：列为 `cell,nGene,nUMI,orig.ident,percent.mito,State`，无 `sample/phase/perc_ribo`；`State` 已含 C/S1–S5，仅作验证 | 2026-10-18 实测 GEO metadata 表头与 State 计数 | D3/D22；GSE206125 metadata | 03/04 | — |
| 2026-10-18 | 单细胞聚类改为**按样本分别聚类**（每线独立 Louvain），与原文 Fig. 1d 的 `C-0…T320-2` 口径一致 | 原文 Methods「individual analysis of the Kuramochi datasets」；审核发现原方案把所有线混聚 | 原文 Methods；Fig. 1d | 03/04 | — |
| 2026-10-18 | 更正归一化口径：原文 TPM 与 Seurat `LogNormalize`（CP10k）不等价；主口径用 CP10k 并声明偏离 | 审核指出"等价"表述错误 | D21；原文 Methods | 03/04/05 | — |
| 2026-10-18 | 模块签名阈值改为 `avg_log2FC > log2(1.8)`（约 0.85），并按 ≥3 signatures + Jaccard 距离 0.25 切树 | 原文"1.8 average fold change"指倍数而非 log2FC；审核发现 `>1.8` 实为 3.5 倍 | D5/D16；原文 Methods | 04 | — |
| 2026-10-18 | CRISPR 设计矩阵按实测 10 个样本列重写：`Samples` 首列、`baseline` 全 1、条件列互斥；T40 仅作敏感性 | 2026-10-18 实测 counts 表头；审核发现原设计条件列全 0 | D10；GSE247686 counts | 08 | — |

> 追加规则：每次 S6 执行、契约变更、参数偏离、数据缺口处理，均在此表追加一行；与 `10_复现质控清单、排错与批判性评述.md` §八 D 表保持同步。
