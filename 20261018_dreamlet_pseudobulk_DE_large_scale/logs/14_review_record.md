# S4 独立审核记录（2026-10-18）

> 依据：`../../REPRO_SPEC.md` L0 §2「S3/S4 必须分离」与 §8.3「四查：事实/代码/引用/结构」。
> 审核方与写作方**不是同一上下文**：写作由主会话完成，审核由两个独立子代理完成（只读，不修改文件）。
> 本记录只登记**发现 → 修正**的对应关系；未修正项明确列出原因。

## 一、审核通道

| 通道 | 职责 | 覆盖 |
|---|---|---|
| `code-reviewer`（Opus） | 统计公式、R API、代码可执行性、内部一致性、诚实性 | 16 篇文档 + `PROJECT_SPEC.md` + PDF 12 页 + dreamlet/variancePartition 源码 |
| `verifier`（Opus） | 数值与事实核对、作者代码参数、补充材料、交叉引用 | 同上 + Supplementary Information（MOESM1，30 页）+ 作者 `dreamlet_analysis` 仓库 |

## 二、发现与修正

### CRITICAL（均已修正）

| # | 发现 | 修正 | 位置 |
|---|---|---|---|
| C1 | 基因 ID 未落地：SYMBOL 查询打在 Ensembl 行名上，PTPRG/MT-/MALAT1 会静默取空 | 02 号聚合后统一转 SYMBOL，落盘 `id_map_<dataset>.tsv`；11 号 QC 前先转；关键基因加 `stopifnot` | 02 §4.2、11 §4.3、PROJECT_SPEC §4 |
| C2 | 残差自由度写成 `n−tr(H)`；EB 误用 `fit$df.residual` | 改为 `rdf = n−2tr(H)+tr(H²)`（`rdf.merMod`）；EB 用 `fit$rdf`；D8 同步 | 04 §2.4/§4、15 §3.4、D8 |
| C3 | `zenith_gsa(coef=)` 参数名错 | 改为 `coefs=` | 07 §4.2/§4.4 |
| C4 | `muscat::aggregateData(cluster_id=, sample_id=)` 参数名错 | 改为 `by=c(...)` | 08 §4.2 |
| C5 | 09/10 的 study-wide BH 跨了多个 coef | 改为 `group_by(coef)` 后再 `p.adjust` | 09、10、12 |
| C6 | 置换在样本级而非受试者级 | 改为 `SubID` 级、pool 内置换后映射回样本 | 08 §4.5 |
| C7 | `sce_psychad_qc/annot.rds` 无生产者 | 11 号补 `saveRDS`；契约表登记 | 11 §4.3/§4.5、PROJECT_SPEC §4 |
| C8 | 使用了不存在的 Hallmark 集合 | p38/synapse 改用 GO（`GOBP_P38MAPK_CASCADE`、`GOBP_SYNAPSE_ASSEMBLY`、`GO1900744`）；Hallmark 仅保留 TNFA | 07、10、12、D9 |

### HIGH（均已修正）

| # | 发现 | 修正 |
|---|---|---|
| H1 | δ 法代数错误（`Var(y)=1/((μ+a)(ln2)²)`） | 改为 `μ/((μ+a)²(ln2)²)`，并注明 `μ≫a` 近似 |
| H2 | 03 号“源码逐行”与实际对象不符 | 标注为摘录/伪代码，保留真实变量名 |
| H3/H4 | 诊断脚本缺输入 | 补 `readRDS`/参数解析/运行前置说明 |
| H5 | 08 的 naive/glmer/绘图是占位或未定义变量 | 补 Python naive、R glmer、CSV 读取；删除 `stop("见 ...")` |
| H6 | 01 安装清单缺包 | 补 `DropletUtils/scran/scDblFinder/bluster/harmony/rtracklayer/Rsamtools/GenomicRanges/fANCOVA/ggrepel/ggpointdensity` |
| H7 | R 版本 4.4/4.5/4.6 不一致 | 统一 R 4.6（Bioc 3.23 要求） |
| H8 | 内存门槛错误地以 naive 为参照 | 改为“相对 Matrix.utils/muscat/pegasus 的内存 <10%”；naive 只比时间 |
| H9 | `pb`/`weights` 保存为 list，契约写 SCE/matrix | 契约改为 list 结构 |
| H10 | `vp_gc` 列名/粒度/脚本名不一致；写入 `raw/` | 改为 per-assay summary；GC 表移到 `processed/benchmark/gc/`；生产者 `06_05` |
| H11 | 方差分数分母漏固定效应；“行内和可 >1” | 按 `calcVarPart` 改为统一归一化（行内和=1） |
| H12 | “每受试者核数”实际按样本算 | 先按 `SubID` 聚合再取均值 |
| H13 | Harmony 成分列名 `V1..V5` 可能取空 | 显式重命名并断言 |
| H14 | 11 号 MT/MALAT1/XIST/UTY 假定 SYMBOL | 先做 ID 转换并断言 |
| H15 | “6 samples/pool”与“6 donors/pool”混用；pool 数无来源 | 统一 6 donors/pool；pool 数改为实测（D18） |
| H16 | `priorWeightsAsCounts=TRUE` 被说成与手工 `Y+0.5` 等价 | 明确“不等价”，仅作敏感性备选 |
| H17 | ŝ 仍被描述为“残差标准差” | 统一为“√(残差标准差) 的预测值” |
| H18 | `plotPercentBars` 传 `DFrame` 会分派失败 | 改为直接传完整 `vpDF` + `genes=/assays=` |

### MEDIUM / LOW（选择性修正）

- M1 VIF 推导补“文库大小相同/可忽略”的前提。
- M2 `theory_sims.Rmd` 出处改为作者仓库，不写“原文补充图”。
- M6 macOS 用 `gtime -v`；脚本按 `uname` 选择。
- M7 Synapse 改用 `synapseutils.syncFromSynapse` 递归下载；去除尖括号占位。
- M9 稀疏内存算术改为 8 B/非零项 ≈58 GB（与原文一致）。
- M13 camera 年份改为 2012。
- M16 pseudobulk 方法数改为 5。
- M17 `delta` 定义补 `G/m2` 缩放。
- M18 注释步骤补 `scoreMarkers` 可执行代码与 `saveRDS`。
- M22 原文 “DYPD” 与 HGNC “DPYD” 的差异在文中说明。
- M23 Gene vs GeneFull 登记为 D19。
- M24 `weights` → `weightsList`。
- L1 内存粗估函数修正单位为 MB。
- L7 校验和改为 `md5 -r`/`md5sum` 跨平台。
- L9 图文件名统一（Fig2A/2B、Fig3G、Fig6A）。
- L10 scDblFinder 与 Scrublet 的算法差异已声明。
- `tools/repro_check.sh` 加执行权限（mode 755）。
- `tools/repro_onboard.sh` 修复 macOS `date -d` 与中文标点紧邻变量的解析崩溃。

## 三、未修正 / 待 S6 实测裁决

| # | 事项 | 原因 | 处理 |
|---|---|---|---|
| R1 | 两级权重到底走 A（显式 `weightsList`）还是 B（`weightsList=NULL`） | 原文 Supplement Eq.16 说有 Poisson 初始化；作者公开代码未显式传；dreamlet 1.10.0 默认全 1 | **D7 规定 A/B 两条都跑**，以能复现 headline 数字者为主路径；S6 填 gap 表 |
| R2 | 作者代码 commit / 补充材料版本 | 已加 `code/01_07_fetch_supplement.sh` 归档步骤，但本次未实际下载 | S6 首次运行时记录 md5 与 commit |
| R3 | `poolID` 的实际数量 | 原文未报告；586/2/6 估算与“~60”不符 | D18：从 public release 实测 |
| R4 | AD 数据 release 与论文 586/299/1.4M 的对应 | public release 滚动更新 | D12：实测后填 gap 表 |
| R5 | Mathys 2023 外部汇总表 | 获取方式未在正文给出 | D20：不可得则登记 gap，不伪造 |
| R6 | 计算基准的绝对耗时/内存 | 硬件不同（512 GB 集群 vs 本机） | D13：只报相对倍数 |

## 四、审核结论

- **事实层**：headline 数字（674K/110、1824、12,472、11.4%、3250、157、525/1899、586/299/1.4M、38.4%、101、932、1037、1.59/2.94e-28、13/14、31 min/26 GB、45 CPU-min）均可在正文定位；作者代码参数（3/3、10/10、5/10、三套公式、三个 contrast）已核对。
- **代码层**：全部 R/Python/bash 代码块通过语法解析；已知不可执行处已改为可运行或明确标注前提。
- **结构层**：`repro_check.sh` 8 项通过；D1–D20 在 L1 索引与 14 号权威表一一对应；表格竖线已转义。
- **结论**：文档集达到“可执行复现方案”的基线；R1（权重路径）与 R2–R6 必须在 S6 用实测数据裁决，已写入 gap 表模板与 D 表。

> **诚实性声明**：本轮审核修正了 8 个 CRITICAL、18 个 HIGH；仍有 6 项依赖 S6 实测。不得把“审核通过”理解为“数值已复现”。
