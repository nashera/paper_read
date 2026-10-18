# S4 独立审核记录

> 项目：`20261018_CRC_myeloid_scRNA_CSF1R_CD40`
> 审核对象：`00`–`15`、`PROJECT_SPEC.md`、`README.md`
> 上位规范：`../REPRO_SPEC.md` v2.1（L0）
> 审核方：独立 `code-reviewer` 通道（与写作方分离）
> 审核日期：2026-10-18
> 审核方式：事实核查（逐条对照原文全文 `/tmp/cell2020_crc_fulltext.txt`）+ 代码/接口核查 + 契约一致性 + 结构检查；`tools/repro_check.sh` 通过。

---

## 一、审核结论

**结论：S4 通过（有条件）。**

- 事实层：关键数值、簇标签、基因列表、LR 对、生存分组与 HR 均与原文一致；M2 基因数已由 26 修正为 **27**（含 `CCL13`）。
- 结构层：16 篇文档头部四要素齐全、①–⑧ 骨架完整、代码围栏配对、无断链、无占位符。
- 接口层：首轮审核发现 8 条严重接口问题（对象类型、平台混用、分母、矩阵方向等），已全部修复并复审闭合。
- 条件：文档中的代码是 **S3 方案级完整代码**，尚未在 S6 实跑；S6 执行前必须按 `14` 号 §六 QC 表与 D 表逐项验证；湿实验按 D19 不执行。

---

## 二、首轮审核的严重问题与修复

| # | 问题 | 修复 | 状态 |
|---|---|---|---|
| 1 | 03 的 `human_annot.rds` 存 Seurat，契约/下游按 SCE 读取 | `03_02b_finalize_annot.R` 统一 `as.SingleCellExperiment()`，补 `logtpm`/`tpm` assay；输出 10x 与 Smart-seq2 两个 SCE | ✅ 闭合 |
| 2 | 08 的 `mouse_annot.rds` 存 Seurat、缺 `cluster_id`/`cell_type`/`mouse_id`/`study`；09–11 接口断裂 | `mouse_annot.rds` 保持 Seurat 粗分；`mouse_myeloid_integrated.rds` 改为带 mM/metadata/logcounts 的 SCE；新增 `08_02b` 生产 `mouse_t_annot.rds`（mT01–mT12）；09/10 读髓系 SCE，11 读 T SCE | ✅ 闭合 |
| 3 | 06/07 用 Smart-seq2 矩阵配 10x 注释对象，`identical(colnames())` 必失败 | 03 新增 `human_smartseq2_annot.rds`；06/07 全部改读该对象；增加细胞交集对齐 | ✅ 闭合 |
| 4 | 04 的 DEG/GSVA 实际落在 10x，却声称 Smart-seq2 | 04_01 新增 `myeloid_human_smartseq2.rds`；markers/DEG/GSVA/M1-M2 改读 Smart-seq2；加平台硬断言 | ✅ 闭合 |
| 5 | 02 把 5,126 既往 T 细胞当“额外合并”，会重复计数 | 明确 10,468 已含 5,126；合并脚本加“已达 10,468 或已有 prior_T_cell 则跳过”的 guard | ✅ 闭合 |
| 6 | 12 的 bulk 矩阵方向写反，signature 全部变 NA | `read_expr()` 统一返回基因 × 样本；`zscore_by_gene()` 按基因跨样本 z-score；`score_sig()` 用 `rownames`；12_03 同步 | ✅ 闭合 |
| 7 | 11 需要 `mT01`–`mT12`，但无生产者 | 08 新增 `08_02b_tcell_annot.R`，输出 `mouse_t_annot.rds`；11 改读该对象 | ✅ 闭合 |
| 8 | Ro/e 文件名/列名与 06 不一致，且缺 `platform` | 03_04 写 `results/ro_e.tsv`，列含 `cluster_id,tissue,observed,expected,ro_e,platform`；06 按 platform 过滤 | ✅ 闭合 |

## 三、主要中等问题与修复

| # | 问题 | 修复 | 状态 |
|---|---|---|---|
| 9 | `hE01_Endothelium-ACKR1` 的 `major_type` 误写 Epithelial | 改为 `Endothelial` | ✅ |
| 10 | CNV correlation 用全局恶性平均，原文要求 same sample | 按 `patient_id` 分别计算恶性参考与 correlation | ✅ |
| 11 | 07 弦图把 SDC2–MMP2 方向画反 | 改为 `MMP2->SDC2`（CAF/内皮 → TAM），输出文件名同步 | ✅ |
| 12 | 11 transition 未按 treatment/timepoint/mouse 分层，列名不符契约 | 按 mouse × treatment × timepoint 分层；输出 `from,to,transition_index,p,n_mice`；Fig 6J 按小鼠展示 | ✅ |
| 13 | 09/10/11 对无 `logcounts` 的人类对象调用 `logcounts()` | 09/10 增加 `get_expr()`（logtpm→logcounts→data）；08 为小鼠 SCE 补 `logcounts` | ✅ |
| 14 | 12_03 `get_sig_safe()` 用 `colnames(scores)` 找基因 | 改为 `rownames(z)`；同步修正矩阵方向 | ✅ |
| 15 | 05 PAGA 读的 h5ad 无生产者；`float(sparse)` 报错 | 05_01 增加 `zellkonverter::writeH5AD()`；PAGA 用 `toarray()` | ✅ |
| 16 | 08 Seurat 5 `merge()` 后未 `JoinLayers()` | 08_01 增加 `JoinLayers()` | ✅ |
| 17 | 04 与 09 都生产 Fig. 4H | 04_05 为唯一生产者；09 只写 `proliferation_scores.tsv` 供核对 | ✅ |
| 18 | 多个输出列与契约不一致 | 逐项对齐：`celltype_signatures`、`correlation_gene_cluster`、`lr_pairs_significant`、`csf1r_de`、`startrac_*`、`th1_cdc1_correlation`、`cdc1_msi` | ✅ |
| 19 | TCR dominant 未按 sample 分组，barcode 会串样 | `group_by(sample_id, barcode, chain)`；`cell_id = sample_id_barcode` | ✅ |
| 20 | 38/36 与 13/4/18/5 被并列，口径不清 | 明确：10x 38 = 13+3+17+5；Smart-seq2 36 = 9+4+18+5；13/4/18/5 为合并标签集合 | ✅ |

## 四、低风险修正

- 07 经验 p 改为 `(b+1)/(B+1)`。
- 13 的 400 mg→400 µg 依据改为 `14 §10.2 第 15 条`。
- 14 Q46 计数容差改 ±5%，Q48 HR 容差改 0.477–0.583（±10%）。
- 00 G1 增加 5–10% 警告区。
- 03 的 `downsample_clusters()` 取消有放回过采样。
- 15 的 Th1-like 不再猜 `hT08?/hT09?`，改为“文献未报告具体簇，用 signature”。
- 12 的 MSI 比例标注为“背景知识（非原文）”。

## 五、复审残留（非阻断）

- 08_01 已把 `study` 规范化为 `Study1`/`Study2` 并校验；若上游 metadata 仍不合法会 `stop`，不会静默。
- 09/10 的“关键行解释”已改为“分母来自 `mouse_cd45_counts.tsv`”。
- 文档代码为 S3 方案级；S6 实跑时必须记录 `sessionInfo`/`pip freeze` 与完整命令，并按 `14` 号 QC 表填实测值。

## 六、审核签署

- 审核方：独立 code-reviewer（S4 通道）
- 结论：**S4 通过（有条件）**；准予进入 S6 执行准备。
- 下一次审核：S6 每模块完成后按 `14` 号 §六 QC 表与 §九 gap 表逐项复核；S6 结题前做一次全项目复审。
