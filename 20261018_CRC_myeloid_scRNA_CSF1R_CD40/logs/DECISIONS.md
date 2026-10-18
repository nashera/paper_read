# 项目决策流水

> 依据：`../REPRO_SPEC.md`（L0 通用规范 v2.1 §10.2）与 `PROJECT_SPEC.md`（L1 项目规范）
> **D 表（D1, D2, …）的权威位置**见 `PROJECT_SPEC.md` §8；本文件只记录**运行时的实际决策与偏离**，不重复 D 表内容。
> 规则：每次改动追加一行，**只增不删**；影响产物的决策必须同步更新 `PROJECT_SPEC.md` 的别名清理表与产物契约表（L0 §3.3 / §4.2）。

| 日期 | 决策 | 触发原因 | 依据（D 编号 / 文献 / 实测） | 影响范围 | 决策人 |
|---|---|---|---|---|---|
| 2026-10-18 | 采用 L0 通用规范 `../REPRO_SPEC.md`（v2.1）作为本项目上位规范，并实例化 L1 `PROJECT_SPEC.md` v1.0 | 项目立项、规范统一 | L0 §0.2、附录 A | 全局 | — |
| 2026-10-18 | 建立目录骨架 `raw/{human_crc_smartseq2,human_crc_10x,tcga_coadread,gtex_colon,mouse_renca,mouse_mc38}`、`processed/`、`results/`、`figures/`、`code/`、`logs/` | 开工前置 | L0 §3.1、附录 B 第 4 步 | 全局 | — |
| 2026-10-18 | 原文 PDF 由 `cell.scRNA-seq和bulk RNA-seq数据推断细胞类群互作网络.pdf` 重命名为 `2020_Zhang_L_Cell_Single-Cell_Analyses_Inform_Mechanisms_of_Myeloid-Targeted_Therapies_in_Colon_Cancer.pdf` 并移入项目目录 | PDF 归档、可读性优先 | L0 §1.2、§3.2；`tools/paper_renamer_cli.py` | 项目根 | — |
| 2026-10-18 | 项目 slug 定为 `CRC_myeloid_scRNA_CSF1R_CD40`；模块编号 00–15（14=QC/D 表，15=精读/背景） | 目录与编号规范 | L0 §1.2–§1.3 | 全局 | — |
| 2026-10-18 | 完成 S3 方案文档 00–15（16 篇）；S4 独立审核进行中 | 方案定稿前不跑数据 | L0 §2（S3/S4/S6） | 全局 | — |
| 2026-10-18 | `results/ro_e.tsv` 统一为 03 号生产、06 号消费；删除模块前缀 `03_`，按共享产物命名 | 06 号契约需要；L0 §3.2/§4.2 | D11；L0 §4.1 | 03, 06 | — |
| 2026-10-18 | 原文 M2 基因列表确认为 **27 个**（含 `CCL13`）；`PROJECT_SPEC.md` §5 由 26 修正为 27 | 04 号逐字核对原文 Methods | 原文 Methods「Define M1/M2」 | 04, 15 | — |
| 2026-10-18 | TCGA/GTEx 生存图（Fig. 4I、5J、7J）与 `survival_cox.tsv` 统一由 **12 号**生产；09/10/11 只提供 cluster 定义或 signature 基因表 | 生产者唯一；避免 09/10 与 12 重复写同一图/表 | L0 §4.2；D17 | 09, 10, 11, 12 | — |
| 2026-10-18 | `activated_dc_signature_genes.tsv` 由 10 号导出，12 号读取；12 号不再重复定义该基因表 | 同上 | L0 §4.2；D17 | 10, 12 | — |
| 2026-10-18 | 修正跨文档章节引用：14 号 QC=§六、gap=§九、批判性评述=§十、D 表=§十一、结题清单=§十二 | 原文引用错位（把 §五 当 QC、§八 当 D 表） | L0 §8.3-4 | 全局 | — |
| 2026-10-18 | 03 号 CNV 代码修正 TxDb Entrez→HGNC symbol 映射；06 号网络代码修正组织特异 `cluster_sets`、无向度数与 `ro_e` 列名 | 代码审查发现的静默失败/口径错误 | L0 §12、附录 C「静默失败」 | 03, 06 | — |
| 2026-10-18 | S4 独立审核（code-reviewer）结论：不通过（8 条严重）；已按 1–8 条修复：03/08 对象统一为 SCE、新增 Smart-seq2 注释对象、04 改用 Smart-seq2 髓系子集、02 修正 5,126 重复计数、12 修正矩阵方向、08 新增 mT01–mT12 生产者、03 Ro/e 输出 `results/ro_e.tsv`（含 platform） | 审核发现人类/小鼠主线接口断裂与静默错误 | L0 §2 S4；`logs/14_review_record.md` 已写，S4 通过（有条件） | 全局 | — |
| 2026-10-18 | 新增契约产物：`processed/human_smartseq2_annot.rds`（03）、`processed/myeloid_human_smartseq2.rds`（04）、`processed/mouse_t_annot.rds`（08）、`processed/mouse_myeloid_integrated.rds` 改为 SCE（08） | 下游按平台/细胞类型取对象，避免混用 10x 与 Smart-seq2 | L0 §4.2 | 03,04,06,07,08,09,10,11 | — |

---
（以下为运行时追加区）
