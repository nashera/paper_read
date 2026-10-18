# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 这个仓库是什么

`paper_read/` 是**文献精读与复现仓库**：以单篇文献为一个复现项目（一个目录），产出可追溯、可核验的复现方案文档（Markdown）与配套代码。**没有传统意义上的 build / test / lint**——"正确性"由两份规范文件 + 自检脚本 + QC 门槛表保证。

所有工作**必须先读** `REPRO_SPEC.md`（根目录，通用规范 v2.1，层级 L0，全仓库强制）。本文件只是快速入口。

## 规范层级（效力：L0 > L1 > L2 > L3）

| 层级 | 文件 | 位置 | 作用 |
|---|---|---|---|
| L0 | `REPRO_SPEC.md` | 仓库根 | 框架与红线：目录/命名/契约/环境/随机性/统计/写作/诚实性，全仓库强制 |
| L1 | `<项目>/PROJECT_SPEC.md` | 项目目录 | L0 在本项目的实例化：文献信息、**产物契约表**、**D 表**（决策登记）、**术语表**、别名清理表 |
| L2 | 项目 QC 门槛表 | 项目 `00`/`NN-1` 号文档 | 各模块 Go/No-Go 数值与允许波动 |
| L3 | 模块文档 `NN_*.md` | 项目目录 | 单模块执行细节 |

**三条总则**：① 诚实性优先于一致性（结果不符如实记录，禁止反向调参凑数）；② 口径唯一（同一概念在一个项目内只有一个定义）；③ 可追溯（任何偏离规范/文献的操作记入 `logs/DECISIONS.md`）。

## 项目结构

```
<YYYYMMDD>_<slug>/          # 一个复现项目 = 一篇文献 = 一个目录
├── PROJECT_SPEC.md         # L1 项目规范
├── <文献>.pdf              # 原文（用 tools/paper_renamer_cli.py 重命名）
├── 00_…NN_*.md             # 模块文档，编号见下
├── raw/ processed/ results/ figures/ code/ logs/
│     # raw=只读原始数据; processed=中间对象; results=结果表; figures=图(pdf+png); code=脚本; logs=日志/DECISIONS.md/校验和
```

### 文档编号约定（固定，不插号）

- `00` 总览与目录（数据资产、全局门槛、流程图）
- `01` 环境与数据获取
- `02 … NN-2` 模块文档（按分析流程：单细胞→调控→轨迹→空间→通讯→bulk/生存→机器学习→延伸→实验验证）
- `NN-1` 复现质控清单、排错与批判性评述（QC 门槛表、gap 表、**D 表权威位置**）
- `NN` 背景知识与方法学原理速查

### 每个模块文档的硬性结构

- **头部四要素**：`目标文献` / `模块对应原文` / `前置依赖` / `复现难度·工时·计算资源`
- **正文骨架 ①–⑧**：① 科学问题与逻辑位置 → ② 背景原理 → ③ 原文逐句拆解 → ④ 逐步可运行代码 → ⑤ 参数解释与调参空间 → ⑥ 预期结果与 QC 检查点 → ⑦ 常见陷阱与排错表 → ⑧ 与文献图表的对应关系表

## 脚本分层（硬性分工）

- **`tools/`**：仓库级**可复用**脚本，逻辑与具体文献无关，不得硬编码项目路径/队列名。写新脚本前先 `ls tools/`。
- **`<项目>/code/`**：项目专用脚本，含项目常量，编号前缀与文档一致（`02_01_read_qc.R`）。
- 项目内脚本被**第二个项目**需要时，剥离项目常量上移到 `tools/`，并在 `REPRO_SPEC.md` 附录 D 登记。

## 常用命令

```bash
# 一致性自检（8 项：结构/围栏/交叉引用/占位符/非法文件名/头部四要素/旧名残留/缺失项标注）
tools/repro_check.sh <项目目录>        # 检查单个项目
tools/repro_check.sh --all             # 扫描仓库下所有含 PROJECT_SPEC.md 的项目
tools/repro_check.sh .                 # 在项目目录内执行

# 新项目脚手架（按 REPRO_SPEC 附录 A+B 生成目录骨架与 PROJECT_SPEC.md，幂等不覆盖）
tools/repro_onboard.sh --dir <项目目录> --paper "<完整文献引用>" [--dry-run]

# 论文 PDF 智能重命名（元数据/首页文本/BibTeX/Crossref/arXiv 多源）
python3 tools/paper_renamer_cli.py <pdf> -t "{year}_{author}_{title}" -n

# 空间分析三件套（长度参数一律 μm，必须显式 --pixel-size-um，缺参直接报错）
python3 tools/spatial_cn_kmeans.py --cells <细胞表> --group-col sample --window-mode radius --radius-um 70 --pixel-size-um 0.5 --k 15 --out-prefix <前缀>
python3 tools/spatial_interaction_test.py --cells <细胞表> --group-col fov --max-dist 30 --dist-unit um --pixel-size-um 0.5 --n-perm 100 --out-prefix <前缀>
python3 tools/patch_proximity.py --cells <细胞表> --patch-by cn --patch-value cn_1 --min-cluster-size 50 --radii-um 5,10,15,20,25 --pixel-size-um 0.5 --out-prefix <前缀>
# 每个脚本自带 --selftest 自检；退出码 0=成功 1=异常 2=校验失败 3=质量断言失败 4=selftest失败 5=缺依赖

# Python 依赖
pip install -r tools/requirements.txt

# 同步（提交 + 推送到 codeup/origin 两个远端；仅在有变更时提交）
./sync.sh               # macOS/Linux
# Windows: sync.ps1（配合 auto.lua + bit.exe + expect.exe）
```

## 关键机制（复现项目的核心纪律）

- **产物契约（L0 §4）**：每个产物**只有一个生产者**，其他模块只能消费。契约表在 `PROJECT_SPEC.md`，字段含 `产物路径 | 生产者(唯一) | 格式/结构(列名·尺度·单位·坐标系) | 消费方 | 生成脚本`。改格式/名 = 改契约表 + 通知全部消费方 + 记 `DECISIONS.md`。
- **D 表（决策登记，L0 §10）**：文献未写清/矛盾处由复现者显式裁决，编号 `D1, D2…` 只增不改号；会改变结论口径的加 ★。权威位置在项目 `NN-1` 号文档（`PROJECT_SPEC.md` 只维护索引）。
- **术语表（唯一口径）**：细胞类型标签、基因/蛋白命名、分组阈值、"阳性 (`<X>`⁺)"判定规则，全项目一个定义。
- **命名规则（L0 §3.2）**：产物/代码/数据文件名 snake_case、禁空格与中文、禁 `< > { } : * ? |`（冒号改连字符）；占位符统一 `<...>` 且落盘必须展开；图表必须 `.pdf`+`.png` 双份。
- **随机种子**：全项目统一常量（建议项目启动日如 `20260101`），写在各脚本首行并在契约表声明。
- **诚实性红线（L0 §12，不可偏离）**：跑不出文献的数不叫失败，改参数直到跑出才叫失败；"没复现出来"≠"文献不可复现"（需 ≥2 条独立路径排除前者）；训练集性能 ≠ 泛化性能；相关性 ≠ 因果；阴性结果必须保留。

## 工作流（六阶段，REPRO_SPEC §2）

S1 立项 → S2 建 `PROJECT_SPEC.md` → S3 写模块文档 → S4 **独立审核**（写作方/审核方必须分离，用 `code-reviewer`/`verifier` 独立通道；四查：事实/代码/引用/结构）→ S5 跑自检 → S6 执行与登记（跑代码→填 QC 表→记 `DECISIONS.md`→出 gap 表）。**S6 才允许跑数据**，方案定稿前不跑分析。

## 已存在的项目

- `20261015_multi-omics_and_machine_learning/` — MASLD→MASH 多组学 + 机器学习
- `20261016_GBM_lactate_NFATc4_SPP1_macrophage/` — 胶质母细胞瘤乳酸/NFATc4–SPP1 轴（文档最完整，可作范例）
- `20261016_GLASS_IDH_glioma_divergent_evolution/` — IDH 胶质瘤分化演化（Nature 2026）
- `20261016_mucosal_melanoma_spatial_multiomics_ICI/` — 黏膜黑色素瘤空间蛋白+空间代谢组学
- `20261018_resistance_continuum_cancer_therapy/` — 肿瘤治疗耐药连续谱（Nature 2024；scRNA/ATAC/WES/CRISPR/PDX）
