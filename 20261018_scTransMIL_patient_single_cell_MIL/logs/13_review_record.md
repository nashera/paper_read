# S4 独立审核记录（2026-10-18）

> 依据：`../../REPRO_SPEC.md` L0 §2「S3/S4 必须分离」与 §8.3「四查：事实/代码/引用/结构」。
> 写作由主会话完成；审核由两个独立只读子代理完成（不修改文件）。
> 本记录只登记**发现 → 修正**的对应关系；未修正项明确列出原因。

## 一、审核通道

| 通道 | 职责 | 覆盖 |
|---|---|---|
| `code-reviewer`（Opus） | 算法/代码正确性、API、内部一致性、诚实性 | 15 篇文档 + PDF 17 页 + 作者仓库 `TencentAILabHealthcare/scTransMIL` + Geneformer/scIB/scGPT/smooth-topk 源码 |
| `verifier`（Opus） | 数值与参数核对、作者代码参数、契约/交叉引用 | 同上 + `1-process_v3.py`/`main.py`/`model_clam.py`/`create_splits_seq.py`/`2-bag_eval.py` |

## 二、发现与修正

### CRITICAL（已修正）

| # | 发现 | 修正 | 位置 |
|---|---|---|---|
| C1 | `gene_id.csv` 按列位置赋值，真实表头是 `Gene stable ID,Gene name` → 映射全为 0 | 改为按列名 `rename`，加 `assert token_id.notna().sum()>0` | 02 §4.1 |
| C2 | tokenizer 代码不可运行（未加 `sys.path`、未传 3 个字典路径、`model_input_size/special_token` 未显式传） | 补 `sys.path.insert`，显式传 4 个路径 + `model_input_size=4096`、`special_token=True` | 02 §4.4 |

### HIGH（已修正）

| # | 发现 | 修正 |
|---|---|---|
| H1 | rank 编码漏掉 `/n_counts*1e4` | 02 §2.3、14 §3.3 补文库归一化 |
| H2 | cell score 用 softmax 前 `A_raw`，与作者 `eval_utils` 保存的 softmax 后权重不符 | 新增 D22：`a_i`（softmax 前）/`attention_weight`（softmax 后）/`cell_score`（AUC 用 softmax 后）；05/09/13/14 统一 |
| H3 | splits 文件格式（长表）与作者/契约（宽表）不兼容 | 05 §4.1 输出 `train/val/test` 三列宽表 + `split_assignments` 长表；契约同步 |
| H4 | `<17` 过滤语义错误（回绕后永不触发；文档还在 tissue 级预过滤） | 删除 tissue 级预过滤；D5/02 注明该检查 inert |
| H5 | tokenized 格式四处不一致（dict vs 两个 .pt；4/5 元组解包） | 统一为 dict(`input_ids`,`attention_mask`,`slide_id`)；02 保存脚本；05 改 5 元组解包 |
| H6 | 03 §2.6 "Stage 1/2 共享 bag embedding" 与代码/ D3 矛盾 | 改为"正文共享 / 代码两个独立任务" |
| H7 | `smooth_topk` 导入不存在 | 改为 `from topk.svm import SmoothTop1SVM` |
| H8 | attention 把位置下标当 token ID | 用 `input_ids` 逐位置映射；过滤 pad |
| H9 | `special_token/model_input_size` 口径矛盾 | 01/02 明确作者实际 `True/4096`；CLS 归约的前提写明 |
| H10 | G2/C1/C3 "全面超越"与作者 50% CSV 矛盾 | G2 限定 100%；50% 时按指标分别报告；13 §9.1 加限定 |
| H11 | 甲状腺 GSE184362 在作者泛癌训练清单中 | 新增 D23：主路径训练排除 GSE184362；09 标 in-corpus case study |
| H12 | D13 自相矛盾（冻结 encoder vs 训最后 2 块）；checkpoint 键名不匹配；`--freeze` 不存在 | D13 改为"冻结 embeddings+前 4 层，训最后 2 块+pooler+aggregator"；10 §4.1 加 key remap；删除 `--freeze` |
| H13 | scIB API 调用错误 | 11 §4.1 改为逐指标函数（nmi/ari/silhouette/silhouette_batch/kBET/graph_connectivity） |
| H14 | scGPT API 不存在（`scgpt.tasks.GeneFormerModel`） | 06 §4.3 改为 `scgpt.tasks.embed_data` / `scgpt.model.GeneformerModel` |
| H15 | 参数量误写 ~30 M | 改为 ~10 M（"30M" 指 Genecorpus-30M 细胞数） |
| H16 | 前向显存 <10 GB 无依据 | 改为按 attention 实现分情况（SDPA 2–4 GB / eager ~26 GB） |
| H17 | 29/30 类的第 30 类被当作"缺失癌种" | D16 改为"29 个癌种标签 + Stage-1 normal 桶（Y_hat=29）" |

### MEDIUM / LOW（选择性修正）

- M1 instance 分类器输入原始 `h`（非 `h_gate`）→ 写入 D3 第四处。
- M2 初始化口径：D10 写明 Geneformer 权重初始化 + 分类头随机；正文"randomly initialized"记录为矛盾。
- M3 契约与文档脚本名统一（01_02/01_04/02_01/02_05/09_03/10_03/10_05）。
- M4 bag_index 加 `dataset` 列 + `02_05_merge_index.py`。
- M5 多分类聚合改为 voting（概率和作敏感性）。
- M6 测试集 sample 数改为"全集 10%"，不再写 1,563。
- M7 训练量子采样按 donor + 癌种分层。
- M8 三个消融改为可运行（MaxPoolCLAM 子类、bag_weight=1.0）。
- M9 gating 指标先 KMeans 再算 NMI/ARI/AMI/HOM。
- M10 关键占位块补成可运行代码（AML 亚型/生存、DEA 导出、embedding 导出）。
- M11 AML 亚型/生存代码补输入定义与基因/患者筛选逻辑。
- M12 09 增加 DEA CSV 产物；11 读取时检查存在。
- M13 20k bags 标注为原文报告值、4M/256 为下界。
- M14 甲状腺细胞数改为"原文未报告"。
- M15 `gc95M` 字典也不在当前 HF 仓库 → D10/01 记录。
- M16 `safe()` 同时用于 `case_id`，并在 D5 记录。
- M17 PCA 改为训练集全局拟合。
- M18 scVI 前检查 X 是否为整数 counts。
- M19 定义 `BagDataset` 返回 `(ids, mask, label, slide_id)`；04/05 统一。
- M20 "工程完整"降级为"核心训练/推理可运行；AML/attention/baseline 需补脚本"。
- M21 生成本审核记录。
- M22 代码固定到 Zenodo commit（01 §4.3）。
- M23 SmoothTop1SVM 在 logits 输入下与 Eq.7 一致；引用改为 Berrada 2018 ICLR。
- M24 Eq.7 实现改为 logits + logsumexp 稳定版。
- M25 冻结策略包含 `pooler`。
- M26 AML 微调强调按 donor 划分（QC 表加交集检查）。
- M27 14 号同步 rank 公式与 SVM 说明。
- M28 13 号注明作者 `metrics_plot/data/*.csv` 可核验。
- M29 00 号注意力表述降级为"相关性线索"。
- M30 loader 接口统一。
- L1 回绕复制次数改为"⌈256/n⌉ 次左右"。
- L2 亚型名补 `-like`。
- L3 低级别胶质瘤例子标为推测。
- L4 显存账按 4 heads × 6 层写清。
- L5 把"以 API 为准"的代码块改为可运行或注明前提。
- L7 DECISIONS D 编号更新为 D1–D23。
- L8 01 号硬件句标注作者 README（非正文）。

## 三、未修正 / 待 S6 实测裁决

| # | 事项 | 原因 | 处理 |
|---|---|---|---|
| R1 | `gf-6L-30M-i2048` checkpoint 与 `gc95M` 字典的最终来源 | 当前 HF 仓库未列出；可能需作者 Drive/Zenodo | D10：先尝试 V1 归档；不可得则用 V1-10M + gc30M 并在方法学声明 |
| R2 | 作者代码确切 commit（Zenodo） | 本次用 `main` 对照；未见 commit hash | S6 首次运行时 checkout Zenodo commit 并记录 |
| R3 | attention 行/列归约方式 | 论文/代码均未给，仓库无提取脚本 | D21：假设 CLS→token，做列求和敏感性 |
| R4 | 29 vs 30 癌种的确切差异 | 需要 Supplementary Data 1 核对 | D16：先按 29+normal 桶报告 |
| R5 | 甲状腺 Patient 5 是否在训练集 | 作者 README 清单包含 GSE184362，但是否含该 sample 未核 | D23：主路径排除 GSE184362，另做不排除对照 |
| R6 | AML 35+6 的 samples/bags 歧义 | 原文 p8/p12 用词不一致 | D14/D5：按 sample 再切 bag，记录实际数 |
| R7 | dbGaP bulk AML 访问 | 受控 | D12：未获批则生存分析标"未复现" |
| R8 | Fig. 2a/2c/3d 的精确数值 | 正文未给；作者 CSV 可核验 | 以 `metrics_plot/data/*.csv` 为准，不编造 |

## 四、审核结论

- **事实层**：29 组原文数值全部可定位；未发现编造图值。
- **代码层**：57 个 Python 块语法通过；作者代码参数（256/1024/<17/回绕/29 类/k=10/lr 5e-5/patience 20）已核对。
- **结构层**：`repro_check.sh` 8 项通过；D1–D23 两表一致；表格竖线已转义。
- **结论**：文档集达到"可执行复现方案"的基线；R1–R8 必须在 S6 用实测数据裁决，已写入 gap 表与 D 表。

> **诚实性声明**：本轮修正了 2 个 CRITICAL、17 个 HIGH；仍有 8 项依赖 S6 实测。不得把"审核通过"理解为"数值已复现"。
