#!/usr/bin/env bash
# ============================================================================
# repro_onboard.sh —— 新复现项目一键脚手架（Onboarding Scaffold）
#
# 上位规范：`REPRO_SPEC.md`（文献复现项目统一规范（通用），L0）
#   · 附录 A  项目规范模板（PROJECT_SPEC.md）→ 本脚本据此逐节生成
#   · 附录 B  新项目落地清单（12 项）        → 生成后打印为"下一步清单"
#   · §1.4    脚本分层：可复用脚本放 `tools/`，不得硬编码项目路径/队列名/论文名
#   · §3.1    项目内目录树（强制）：raw/ processed/ results/ figures/ code/ logs/
#   · §6.1    随机种子常量 = 项目启动日 YYYYMMDD
#   · §10.2   `logs/DECISIONS.md` 决策流水模板
#
# 用途：给一个新项目一键生成**符合 L0 的目录骨架与文档骨架**。
#       只生成结构与占位，**绝不编造文献内容**；已存在的文件一律不覆盖。
#
# 用法：
#   tools/repro_onboard.sh --dir <项目目录> --paper "<完整文献引用>" \
#                          [--date YYYY-MM-DD] [--raw-dirs <a,b,...>] \
#                          [--dry-run] [--force]
#   tools/repro_onboard.sh --help
#
#   参数：
#     --dir   <路径>        要创建/初始化的项目目录（相对或绝对）。
#                           若已存在且非空 → 报错退出（退出码 2），除非给 --force。
#     --paper "<引用>"      文献完整引用（写入 PROJECT_SPEC.md 与 00 号总览头部）。
#     --date  YYYY-MM-DD    项目启动日期，默认今天。用于 00 号文档头部与 set.seed 常量提示。
#     --raw-dirs <a,b,...>  raw/ 下的子目录列表（默认只建 raw/ 本身）。禁止空格与 L0 §3.2 禁用字符。
#     --dry-run             只打印将要执行的动作，不落盘。
#     --force               允许在非空目录中**补齐缺失文件**（绝不覆盖已存在的文件）。
#     -h, --help            打印本帮助并退出（退出码 0）。
#
# 示例：
#   # 预演（不落盘）
#   tools/repro_onboard.sh --dir 20261101_foo_gbm --paper "Doe J, et al. *A title.* Nature (2026) 1:2." \
#                          --date 2026-11-01 --raw-dirs geo,spatial --dry-run
#   # 实跑
#   tools/repro_onboard.sh --dir 20261101_foo_gbm --paper "Doe J, et al. *A title.* Nature (2026) 1:2." \
#                          --raw-dirs geo,spatial
#   # 补齐缺失文件（不覆盖已有）
#   tools/repro_onboard.sh --dir 20261101_foo_gbm --paper "..." --force
#
# 生成物：
#   <dir>/raw/{<--raw-dirs 各项>}  processed/  results/  figures/  code/  logs/
#   <dir>/PROJECT_SPEC.md             （按 L0 附录 A §1–§11 逐节生成）
#   <dir>/logs/DECISIONS.md           （按 L0 §10.2 模板，预填 2 行）
#   <dir>/00_目录与复现总览.md          （最小骨架，明确标注待填充）
#
# 退出码：0 = 成功；1 = 未预期异常；2 = 参数/环境校验失败
# 依赖：bash、date、mkdir、cat、sed、tr、printf（无 jq / python 依赖）
# ============================================================================

set -euo pipefail

PROG="$(basename -- "$0")"

E_OK=0
E_UNEXPECTED=1
E_USAGE=2

# ---------------------------------------------------------------------------
# 0. 参数变量
# ---------------------------------------------------------------------------
DIR=""
PAPER=""
PROJECT_DATE=""
RAW_DIRS=""
DRY_RUN=0
FORCE=0
SHOW_HELP=0

RAW_DIRS_NORM=""      # 校验并规范化后的逗号列表
RAW_DESC="raw/"       # 供文档使用的 raw 目录写法
ROOT_ABS=""           # 目标目录绝对路径
SPEC_PATH=""          # L0 规范文件路径
SPEC_REF=".."         # 写入文档的"上位规范"引用写法
SPEC_VERSION="版本见该文件"
SEED=""
WRITE_IT=0            # 由 should_write 设置的写盘开关
WRITTEN=0
SKIPPED=0

# ---------------------------------------------------------------------------
# 1. 输出与退出
# ---------------------------------------------------------------------------
say()  { printf '%s\n' "$*"; }
step() { printf '\n[%s] %s\n' "$1" "$2"; }
die_usage() {
  printf '%s: 参数错误：%s\n' "$PROG" "$*" >&2
  printf '试运行 `%s --help` 查看用法。\n' "$PROG" >&2
  exit "$E_USAGE"
}
die_env() {
  printf '%s: 环境错误：%s\n' "$PROG" "$*" >&2
  exit "$E_USAGE"
}
trap 'rc=$?; printf "%s: 未预期异常（第 %s 行，退出码 %s）\n" "$PROG" "$LINENO" "$rc" >&2; exit '"$E_UNEXPECTED" ERR

usage() {
  cat <<'USAGE'
repro_onboard.sh —— 新复现项目一键脚手架
上位规范：REPRO_SPEC.md（通用规范 L0）附录 A（PROJECT_SPEC 模板）/ 附录 B（落地清单）/ §1.4 / §3.1 / §10.2

用法：
  tools/repro_onboard.sh --dir <项目目录> --paper "<完整文献引用>" \
                         [--date YYYY-MM-DD] [--raw-dirs <a,b,...>] \
                         [--dry-run] [--force]
  tools/repro_onboard.sh --help

参数：
  --dir <路径>          必填。要创建/初始化的项目目录（相对或绝对路径）。
                        若目录已存在且非空 → 报错退出（退出码 2），除非给 --force。
  --paper "<引用>"      必填。文献完整引用，写入 PROJECT_SPEC.md 与 00 号总览的头部。
  --date YYYY-MM-DD     可选。项目启动日期，默认今天。
                        用于 00 号文档头部，并推导全项目随机种子常量（YYYYMMDD，L0 §6.1）。
  --raw-dirs <a,b,...>  可选。raw/ 下的子目录列表（逗号分隔），默认只建 raw/ 本身。
                        名称禁止空格与 L0 §3.2 禁用字符 < > { } : * ? | 及斜杠。
  --dry-run             可选。只打印将要执行的动作，不写任何文件。
  --force               可选。允许在非空目录中**补齐缺失文件**；已存在的文件一律跳过，绝不覆盖。
  -h, --help            打印本帮助并退出。

生成物：
  <dir>/raw/{<--raw-dirs 各项>} processed/ results/ figures/ code/ logs/
  <dir>/PROJECT_SPEC.md       按 L0 附录 A 模板逐节生成（§1–§11）
  <dir>/logs/DECISIONS.md     按 L0 §10.2 模板生成，预填 2 行
  <dir>/00_目录与复现总览.md    最小骨架（头部占位 + 章节标题 + 待填充标记），不编造任何文献内容

示例：
  tools/repro_onboard.sh --dir 20261101_foo_gbm \
                         --paper "Doe J, et al. *A title.* Nature (2026) 1:2. DOI: 10.1000/xyz" \
                         --date 2026-11-01 --raw-dirs geo,spatial --dry-run
  tools/repro_onboard.sh --dir 20261101_foo_gbm --paper "..." --raw-dirs geo,spatial
  tools/repro_onboard.sh --dir 20261101_foo_gbm --paper "..." --force

退出码：
  0  成功
  1  未预期异常
  2  参数/环境校验失败（缺 --dir/--paper、日期非法、目录非空且无 --force 等）

依赖：bash、date、mkdir、cat、sed、tr、printf（无 jq / python 依赖）
USAGE
}

# ---------------------------------------------------------------------------
# 2. 参数解析与校验
# ---------------------------------------------------------------------------
require_value() {
  # $1 = 选项名；$2 = 剩余参数个数（应为 $#）
  if [ "$2" -lt 2 ]; then
    die_usage "选项 $1 缺少参数值"
  fi
}

parse_args() {
  while [ $# -gt 0 ]; do
    case "$1" in
      -h|--help)   SHOW_HELP=1; shift ;;
      --dir)       require_value "$1" "$#"; DIR="$2"; shift 2 ;;
      --dir=*)     DIR="${1#--dir=}"; shift ;;
      --paper)     require_value "$1" "$#"; PAPER="$2"; shift 2 ;;
      --paper=*)   PAPER="${1#--paper=}"; shift ;;
      --date)      require_value "$1" "$#"; PROJECT_DATE="$2"; shift 2 ;;
      --date=*)    PROJECT_DATE="${1#--date=}"; shift ;;
      --raw-dirs)  require_value "$1" "$#"; RAW_DIRS="$2"; shift 2 ;;
      --raw-dirs=*) RAW_DIRS="${1#--raw-dirs=}"; shift ;;
      --dry-run)   DRY_RUN=1; shift ;;
      --force)     FORCE=1; shift ;;
      --)          shift; break ;;
      -*)          die_usage "未知选项：$1" ;;
      *)           die_usage "不接受位置参数：$1（请用 --dir / --paper）" ;;
    esac
  done
  if [ $# -gt 0 ]; then
    die_usage "不接受位置参数：$1（请用 --dir / --paper）"
  fi
}

abs_path() {
  local p="$1" r=""
  if command -v realpath >/dev/null 2>&1; then
    if r="$(realpath -m -- "$p" 2>/dev/null)" && [ -n "$r" ]; then
      printf '%s\n' "$r"; return 0
    fi
  fi
  case "$p" in
    /*) printf '%s\n' "$p" ;;
    *)  printf '%s/%s\n' "$(pwd)" "$p" ;;
  esac
}

validate_args() {
  if [ "$SHOW_HELP" -eq 1 ]; then usage; exit "$E_OK"; fi

  if [ -z "$DIR" ]; then die_usage "缺少必填参数 --dir（用法见 --help）"; fi
  if [ -z "$PAPER" ]; then die_usage "缺少必填参数 --paper（用法见 --help）"; fi

  # --- 目录名合法性（L0 §3.2：禁空格与 < > { } : * ? |） --------------------
  local base
  base="$(basename -- "$DIR")"
  if [ -z "$base" ] || [ "$base" = "/" ] || [ "$base" = "." ] || [ "$base" = ".." ]; then
    die_usage "--dir 不是合法的项目目录名：$DIR"
  fi
  case "$base" in
    *" "*) die_usage "--dir 的目录名不得含空格（L0 §3.2）：$base" ;;
    *[\<\>\{\}\:\*\?\|]*) die_usage "--dir 的目录名含 L0 §3.2 禁用字符：$base" ;;
  esac

  # --- 日期 ----------------------------------------------------------------
  if [ -z "$PROJECT_DATE" ]; then
    if ! PROJECT_DATE="$(date +%Y-%m-%d)"; then
      die_env "无法获取当前日期（date 命令失败），请显式提供 --date YYYY-MM-DD"
    fi
  fi
  case "$PROJECT_DATE" in
    [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]) ;;
    *) die_usage "--date 格式应为 YYYY-MM-DD，收到：$PROJECT_DATE" ;;
  esac
  # GNU date 用 -d，BSD/macOS date 用 -j -f（L0 §1.4：tools/ 脚本需跨平台可运行）
  local norm=""
  if date -d "$PROJECT_DATE" +%Y-%m-%d >/dev/null 2>&1; then
    norm="$(date -d "$PROJECT_DATE" +%Y-%m-%d 2>/dev/null)"
  elif date -j -f "%Y-%m-%d" "$PROJECT_DATE" +%Y-%m-%d >/dev/null 2>&1; then
    norm="$(date -j -f "%Y-%m-%d" "$PROJECT_DATE" +%Y-%m-%d 2>/dev/null)"
  else
    die_usage "--date 不是合法日历日期：$PROJECT_DATE"
  fi
  if [ "$norm" != "$PROJECT_DATE" ]; then
    die_usage "--date 不是合法日历日期：$PROJECT_DATE"
  fi
  if ! SEED="$(printf '%s' "$PROJECT_DATE" | tr -d '-')"; then
    die_env "无法由日期推导随机种子常量，请检查 --date"
  fi

  # --- raw 子目录 ----------------------------------------------------------
  RAW_DIRS_NORM=""
  if [ -n "$RAW_DIRS" ]; then
    local stripped="" item="" accum="" _ifs="$IFS"
    # 只容许"逗号两侧"的空白（如 "a, b"），项内空白一律报错——不静默改写用户输入
    if ! stripped="$(printf '%s' "$RAW_DIRS" | sed -e 's/[[:space:]]*,[[:space:]]*/,/g')"; then
      die_env "无法规范化 --raw-dirs"
    fi
    case "$stripped" in
      *[[:space:]]*)
        die_usage "--raw-dirs 的目录名不得含空格（L0 §3.2）：$RAW_DIRS" ;;
    esac
    IFS=','
    for item in $stripped; do
      IFS="$_ifs"
      if [ -z "$item" ]; then
        die_usage "--raw-dirs 含空项（检查多余的逗号）：$RAW_DIRS"
      fi
      case "$item" in
        .|..)  die_usage "--raw-dirs 不允许 . 或 ..：$item" ;;
        -*)    die_usage "--raw-dirs 的项不得以 - 开头：$item" ;;
        */*)   die_usage "--raw-dirs 只接受单层名称（不含 /）：$item" ;;
        *[\<\>\{\}\:\*\?\|]*) die_usage "--raw-dirs 含 L0 §3.2 禁用字符：$item" ;;
        *[[:space:]]*) die_usage "--raw-dirs 的项不得含空格：$item" ;;
      esac
      if [ -z "$accum" ]; then accum="$item"; else accum="$accum,$item"; fi
      IFS=','
    done
    IFS="$_ifs"
    RAW_DIRS_NORM="$accum"
  fi
  if [ -n "$RAW_DIRS_NORM" ]; then
    RAW_DESC="raw/{$RAW_DIRS_NORM}"
  fi
}

resolve_paths() {
  ROOT_ABS="$(abs_path "$DIR")"

  # 仓库根需与 ROOT_ABS 用同一套规范化（realpath 会解析符号链接；若仓库根本身是
  # 软链挂载，两侧不做同一套解析，相对写法判定会失效）
  local script_dir="" repo_root=""
  script_dir="$(abs_path "$(dirname -- "$0")")"
  repo_root="$(dirname -- "$script_dir")"
  SPEC_PATH="$repo_root/REPRO_SPEC.md"

  # 解析 L0 头部 "**版本**：**vX.Y**"；读不到则退化写法，不阻断
  if [ -f "$SPEC_PATH" ]; then
    local v=""
    v="$(sed -n 's/.*\*\*版本\*\*：\*\*\(v[0-9][0-9.]*\)\*\*.*/\1/p' "$SPEC_PATH" | head -1)"
    if [ -n "$v" ]; then SPEC_VERSION="$v"; fi
  fi

  # 项目目录恰为仓库根的一级子目录时用相对写法（L0 附录 A 模板写法），否则用绝对路径
  if [ "$(dirname -- "$ROOT_ABS")" = "$repo_root" ]; then
    SPEC_REF='../REPRO_SPEC.md'
  else
    SPEC_REF="$SPEC_PATH"
  fi

  # --- 目标目录状态检查 ----------------------------------------------------
  if [ -e "$ROOT_ABS" ] && [ ! -d "$ROOT_ABS" ]; then
    die_env "--dir 已存在且不是目录：$ROOT_ABS"
  fi
  if [ -d "$ROOT_ABS" ]; then
    local has_content=0
    if [ -n "$(ls -A -- "$ROOT_ABS" 2>/dev/null || true)" ]; then has_content=1; fi
    if [ "$has_content" -eq 1 ] && [ "$FORCE" -eq 0 ]; then
      die_env "目标目录已存在且非空：$ROOT_ABS
  若只想补齐缺失文件（绝不覆盖已有文件），请加 --force。"
    fi
  fi
}

# ---------------------------------------------------------------------------
# 3. 落盘原语（幂等、不覆盖）
# ---------------------------------------------------------------------------
rel() {
  case "$1" in
    "$ROOT_ABS"/*) printf '%s\n' "${1#"$ROOT_ABS"/}" ;;
    *)             printf '%s\n' "$1" ;;
  esac
}

make_dir() {
  local d="$1" shown=""
  shown="$(rel "$d")"
  if [ -d "$d" ]; then
    printf '  ·  目录已存在：%s\n' "$shown"
  elif [ "$DRY_RUN" -eq 1 ]; then
    printf '  [dry-run] 将创建目录：%s/\n' "$shown"
  else
    mkdir -p -- "$d"
    WRITTEN=$((WRITTEN + 1))
    printf '  📁 创建目录：%s/\n' "$shown"
  fi
}

should_write() {
  local f="$1" shown=""
  shown="$(rel "$f")"
  if [ -e "$f" ]; then
    SKIPPED=$((SKIPPED + 1))
    printf '  ⏭️  已存在，跳过（不覆盖）：%s\n' "$shown"
    WRITE_IT=0
  elif [ "$DRY_RUN" -eq 1 ]; then
    printf '  [dry-run] 将创建文件：%s\n' "$shown"
    WRITE_IT=0
  else
    printf '  ✍️  创建文件：%s\n' "$shown"
    WRITE_IT=1
  fi
}

# ---------------------------------------------------------------------------
# 4. 生成 PROJECT_SPEC.md（L0 附录 A 模板逐节）
# ---------------------------------------------------------------------------
gen_project_spec() {
  local out="$ROOT_ABS/PROJECT_SPEC.md"
  should_write "$out"
  if [ "$WRITE_IT" -ne 1 ]; then return 0; fi

  cat > "$out" <<EOF
# PROJECT_SPEC —— 项目规范（L1）

> **上位规范**：\`$SPEC_REF\`（通用规范 **$SPEC_VERSION**，层级 **L0**，全仓库强制）
> **本文件层级**：**L1**（通用规范在本项目的实例化；不得违背 L0 红线，可填写 L0 未规定的项目特有内容）
> **本项目**：$PAPER
> **项目目录**：\`$ROOT_ABS\`
> **启动日期**：${PROJECT_DATE}　**本文件版本**：**v1.0**

> **效力**：L0 \`$SPEC_REF\` > **本文件 L1** > 项目 QC 门槛表（L2）> 各模块文档（L3）；冲突时以高层级为准（L0 §0.2）。
> **骨架说明**：本文件由 \`tools/repro_onboard.sh\` 依 L0 附录 A 模板生成，**只含结构与占位**。
> 凡标注「（复现者填写）」处必须依 L0 §8 补齐，**不得以空骨架交付**。
> 除上方由命令行参数确定的三项（文献引用、项目目录、启动日期）外，本文件**不含任何对该文献的推论**；
> 凡标「文献报告」的数值必须能在原文中定位（L0 §8.3 红线 1）。

---

## 1. 项目定义

- 一句话目标 / 科学问题：（复现者填写）
- 交付物清单（文档 NN 篇、代码、产物）：（复现者填写）
- 模块划分与编号表（编号规则见 L0 §1.3：\`00\` 总览与目录 / \`01\` 环境与数据获取 / \`02\` 起模块文档 / \`NN-1\` 质控与批判性评述 / \`NN\` 背景知识与方法学速查；**编号只顺延、不插号**）：

| 编号 | 文档标题 | 对应原文位置（Fig./Table） | 状态 |
|---|---|---|---|
| 00 | 目录与复现总览（\`00_目录与复现总览.md\`） | — | 骨架已生成 |
| 01 | 环境与数据获取 | （复现者填写） | 未开始 |
| 02 | （复现者填写，顺延编号） | | 未开始 |

## 2. 目录与命名（本项目的具体化）

- 项目 slug 与目录名：\`$ROOT_ABS\`
- 数据源目录细分：\`$RAW_DESC\`（若此处仅为 \`raw/\`，请按实际数据源补写子目录名）
- 命名特例（如基因集文件、图像文件的特殊约定）：（复现者填写）
- 本项目遵循 L0 §3.2：\`<前缀>_<语义描述>[__<限定>].<扩展名>\`；产物一律 \`snake_case\` 且**不含中文、不留尖括号占位符**；
  图同时出 \`.pdf\` + \`.png\`（dpi = 300）。

## 3. 别名清理表

| 历史变体 | ✅ 规范名 | 出处 |
|---|---|---|
| 无 | 无 | — |

> 规则（L0 §3.3）：改名必须**同时**完成 ① 改所有引用处 ② 更新 §4 契约表 ③ 记 \`logs/DECISIONS.md\`。

## 4. 产物契约表

| 产物路径 | 生产者(唯一) | 格式/结构（列名·尺度·单位·坐标系） | 消费方 | 生成脚本 |
|---|---|---|---|---|
| （复现者填写：核心共享对象 + 各模块产物，逐条登记） | | | | |

> 纪律（L0 §4.2）：生产者唯一；格式即契约（列名、维度、log/非 log、单位、坐标系必须写清）；
> **禁止静默改格式**；跨语言交接用 \`.h5ad\`/\`.csv\`；上游缺失时明确报错而非自动猜路径。

## 5. 术语表

| 术语 | 本项目定义（含阳性判定规则、阈值） |
|---|---|
| （复现者填写） | |

> 任何「阳性 (<X>⁺)」标签必须给出判定规则（表达阈值 or 打分阈值），全项目唯一（L0 §0.3 口径唯一、§8.3-5）。

## 6. 环境与版本钉死

| 组件 | 版本 | 出处 | 偏离说明 |
|---|---|---|---|
| 主语言 / 运行时 | （复现者填写） | 文献指定 or 我方选择 | |
| 关键包 | （复现者填写） | | |
| 参考数据库（注释/motif/基因集…） | （复现者填写版本 + 日期） | | |

> 要求（L0 §5）：版本钉死并提供环境文件（\`conda env export\` / \`renv.lock\` / \`requirements.txt\`）；
> 参考数据库必须记录**版本与日期**；每脚本结束写会话快照并记录**完整执行命令**。

## 7. 统计与随机性

- **随机种子常量**：\`$SEED\`（= 项目启动日 $PROJECT_DATE 的 YYYYMMDD 形式，L0 §6.1）
  - R：\`set.seed($SEED)\`　Python：\`random.seed($SEED)\` / \`np.random.seed($SEED)\`
  - 写在各脚本首行，并在 §4 契约表中声明。
- 必须固定种子的环节（L0 §6.2）：降维与嵌入、聚类、重采样、置换检验、模型、图推断。
- 已知不可完全固定：分布式/并行实现（dask、多线程 BLAS）→ 报告中必须声明。
- 项目特有统计约定：（复现者填写；默认遵循 L0 §7——多重比较 BH(FDR)、生存报 HR + 95% CI + n、
  空间统计置换零模型须保持空间自相关、C-index 与 time-dependent AUC 不可混用）。

## 8. 决策登记表（D 表）

| # | 决策项 | 涉及模块 | 裁决(主路径) | 备选/敏感性 | 依据 | ★ |
|---|---|---|---|---|---|---|
| （复现者填写：自 D1 起顺延，只增不改号） | | | | | | |

> ★ = **改变结论口径**，必须写进论文方法学（L0 §10.1）。
> 每条决策须能追溯到「文献缺失 / 矛盾」的证据（引原文句，或注明原文无该信息，见 L0 §9）。
> **D 表权威位置（二选一，全项目唯一）**：☐ 本文件 §8　☐ 项目质控文档 \`NN-1\`。

## 9. QC 门槛

- **数值层**：文献报告的关键数值 —— 与文献一致（允许平台/抽样差异，**默认 ±10%**，逐项声明）。
- **判据层**：✅ 通过 → 推进；⚠️ 小幅偏差 → 排查实现细节后记录放行；❌ 显著偏差 → **第二套独立实现复核**，仍不符则登记 gap 表。
- **结题层**：gap 表「文献值 / 我方复现值 / 是否一致 / 差异原因 / 处理方式」。
- 全局 Go/No-Go 数值：（复现者填写，并指向项目质控文档的检查点表）
> **禁止**为通过门槛反复调参（L0 §10.3、§12 红线 1）。

## 10. 偏离声明（对 L0 的非红线偏离）

| L0 条目 | 偏离内容 | 理由 | 影响范围 |
|---|---|---|---|
| 无 | 无 | — | — |

> L0 §12 **诚实性红线不可偏离**；任何非红线偏离必须在此登记（条目、理由、影响范围）。

## 11. 修订记录

| 版本 | 日期 | 修改内容 | 修改人 |
|---|---|---|---|
| v1.0 | $PROJECT_DATE | 由 \`tools/repro_onboard.sh\` 依 L0 附录 A 模板生成项目骨架 | \`tools/repro_onboard.sh\` |
EOF
  WRITTEN=$((WRITTEN + 1))
}

# ---------------------------------------------------------------------------
# 5. 生成 logs/DECISIONS.md（L0 §10.2 模板）
# ---------------------------------------------------------------------------
gen_decisions() {
  local out="$ROOT_ABS/logs/DECISIONS.md"
  should_write "$out"
  if [ "$WRITE_IT" -ne 1 ]; then return 0; fi

  cat > "$out" <<EOF
# 项目决策流水

> 依据：\`$SPEC_REF\`（L0 通用规范 $SPEC_VERSION §10.2）与 \`PROJECT_SPEC.md\`（L1 项目规范）
> **D 表（D1, D2, …）的权威位置**见 \`PROJECT_SPEC.md\` §8；本文件只记录**运行时的实际决策与偏离**，不重复 D 表内容。
> 规则：每次改动追加一行，**只增不删**；影响产物的决策必须同步更新 \`PROJECT_SPEC.md\` 的别名清理表与产物契约表（L0 §3.3 / §4.2）。

| 日期 | 决策 | 触发原因 | 依据（D 编号 / 文献 / 实测） | 影响范围 | 决策人 |
|---|---|---|---|---|---|
| ${PROJECT_DATE} | 采用 L0 通用规范 \`${SPEC_REF}\`（${SPEC_VERSION}）作为本项目上位规范，并实例化 L1 \`PROJECT_SPEC.md\` v1.0 | 项目立项、规范统一 | L0 §0.2、附录 A | 全局 | — |
| $PROJECT_DATE | 建立目录骨架 \`$RAW_DESC\`、\`processed/\`、\`results/\`、\`figures/\`、\`code/\`、\`logs/\` | 开工前置 | L0 §3.1、附录 B 第 4 步 | 全局 | — |

---
（以下为运行时追加区）
EOF
  WRITTEN=$((WRITTEN + 1))
}

# ---------------------------------------------------------------------------
# 6. 生成 00_目录与复现总览.md（最小骨架，不编造内容）
# ---------------------------------------------------------------------------
gen_overview() {
  local out="$ROOT_ABS/00_目录与复现总览.md"
  should_write "$out"
  if [ "$WRITE_IT" -ne 1 ]; then return 0; fi

  cat > "$out" <<EOF
# 00 目录与复现总览

> **目标文献**：$PAPER
> **项目目录**：\`$ROOT_ABS\`
> **启动日期**：$PROJECT_DATE
> **上位规范**：\`${SPEC_REF}\`（通用规范 ${SPEC_VERSION}，L0）；本文件属 L3 文档，只可引用 L1 的口径

---

## 0. 头部信息（按 L0 §8.1 四要素）

> **目标文献**：$PAPER
> **模块对应原文**：（复现者填写：本总览索引的 Fig./Table 范围）
> **前置依赖**：（复现者填写；总览为全项目入口，通常填「无」）
> **复现难度 / 预计工时 / 计算资源**：（复现者填写：★★★☆☆ / x 人日 / 内存·CPU·GPU）

## 1. 文档清单

（复现者填写：\`00\`–\`NN\` 的编号、标题、状态；编号规则见 L0 §1.3，编号只顺延不插号）

## 2. 科学问题

（复现者填写：原文的核心科学问题与主结论链）

## 3. 全流程

（复现者填写：模块依赖与执行顺序，可用 mermaid 流程图）

## 4. 图表 → 文档映射

| 原文图/表 | 对应模块文档 | 复现产物 |
|---|---|---|
| （复现者填写） | | |

## 5. 数据资产

| 数据集 / 队列 | 样本数 | 用途 | 获取方式 | 本地路径 |
|---|---|---|---|---|
| （复现者填写） | | | | |

## 6. 算力与工期

（复现者填写：内存 / CPU / GPU / 预计人日）

## 7. 全局门槛

（复现者填写：全局 Go/No-Go；细表见项目 L2 质控文档）

## 8. 工程约定

（复现者填写；默认遵循 L0：\`raw/\` 只读、脚本幂等、图出 \`.pdf\` + \`.png\` 双份、每张图可追溯到「脚本 + 输入产物 + 参数」）

---

<!-- 骨架占位：本文件由 tools/repro_onboard.sh 生成，上面各章节须按 L0 §8 写作规范由复现者补齐。
     正文中「（复现者填写）」即待完成标记；本注释刻意不含 tools/repro_check.sh 第 4 项
     所扫描的那几个关键词，以免刚生成的新骨架被误判为"占位符残留"并返回退出码 1。 -->
EOF
  WRITTEN=$((WRITTEN + 1))
}

# ---------------------------------------------------------------------------
# 7. 目录骨架
# ---------------------------------------------------------------------------
create_skeleton() {
  step 1/5 "目录骨架（L0 §3.1 强制）"
  make_dir "$ROOT_ABS"
  make_dir "$ROOT_ABS/raw"
  if [ -n "$RAW_DIRS_NORM" ]; then
    local _ifs="$IFS" sub=""
    IFS=','
    for sub in $RAW_DIRS_NORM; do
      IFS="$_ifs"
      make_dir "$ROOT_ABS/raw/$sub"
      IFS=','
    done
    IFS="$_ifs"
  fi
  make_dir "$ROOT_ABS/processed"
  make_dir "$ROOT_ABS/results"
  make_dir "$ROOT_ABS/figures"
  make_dir "$ROOT_ABS/code"
  make_dir "$ROOT_ABS/logs"
}

print_next_steps() {
  local state="均未完成"
  if [ "$DRY_RUN" -eq 1 ]; then state="dry-run 未落盘"; fi
  step 5/5 "下一步清单（L0 附录 B，${state}）"
  cat <<'NEXT'
  - [ ]  1. 目录命名 <YYYYMMDD>_<slug>，原文 PDF 归档并去掉非法字符
  - [ ]  2. 读全文，拆模块（01 … NN），确定 00 号骨架
  - [ ]  3. 按附录 A 补全 PROJECT_SPEC.md，先写 §1/§2/§6
  - [ ]  4. 目录骨架 + logs/DECISIONS.md；先 `ls ../tools/` 看有哪些可复用脚本（§1.4：复用优先，禁止把 tools/ 脚本复制进项目）
  - [ ]  5. 数据资产清单（数据集/样本数/用途/获取方式）+ raw/ 校验和
  - [ ]  6. 逐模块写方案文档（头部四要素 + ①–⑧ 骨架）
  - [ ]  7. 补全产物契约表 §4（生产者唯一，格式写细）
  - [ ]  8. 逐条梳理文献缺失/矛盾 → 形成 D 表（★ 标注口径类）
  - [ ]  9. 建 QC 门槛表（数值层 / 判据层 / 结题层 gap 表）
  - [ ] 10. 独立审核（事实/代码/引用/结构四查），修正后跑 tools/repro_check.sh
  - [ ] 11. 开跑：按模块执行，逐模块填 QC 表 + 记 DECISIONS.md
  - [ ] 12. 结题：填 gap 表 + 写批判性评述（证据等级 / 扩展方向）
NEXT
}

# ---------------------------------------------------------------------------
# 8. 主流程
# ---------------------------------------------------------------------------
main() {
  parse_args "$@"
  validate_args
  resolve_paths

  say "repro_onboard.sh —— 新复现项目一键脚手架（L0 附录 A/B 落地）"
  say "上位规范：${SPEC_PATH}（${SPEC_VERSION}）"
  say "────────────────────────────────────────────────────────────"
  say "目标目录：$ROOT_ABS"
  say "文献引用：$PAPER"
  say "启动日期：${PROJECT_DATE}（随机种子常量 ${SEED}）"
  if [ -n "$RAW_DIRS_NORM" ]; then
    say "raw 子目录：$RAW_DIRS_NORM"
  else
    say "raw 子目录：（无，仅建 raw/ 本身）"
  fi
  if [ "$DRY_RUN" -eq 1 ]; then
    say "模式：--dry-run（只打印，不落盘）"
  elif [ "$FORCE" -eq 1 ]; then
    say "模式：--force（非空目录中补齐缺失文件；已存在文件一律跳过）"
  else
    say "模式：常规（不覆盖任何已存在文件）"
  fi
  say "────────────────────────────────────────────────────────────"

  create_skeleton

  step 2/5 "项目规范 PROJECT_SPEC.md（L0 附录 A，§1–§11）"
  gen_project_spec

  step 3/5 "决策流水 logs/DECISIONS.md（L0 §10.2）"
  gen_decisions

  step 4/5 "总览骨架 00_目录与复现总览.md（L0 §1.3）"
  gen_overview

  print_next_steps

  printf '\n────────────────────────────────────────────────────────────\n'
  if [ "$DRY_RUN" -eq 1 ]; then
    printf '结果：✅ dry-run 完成（未写入任何文件）\n'
  else
    printf '结果：✅ 完成（新建 %s 项，跳过已存在 %s 项）\n' "$WRITTEN" "$SKIPPED"
  fi
  printf '建议下一步：bash tools/repro_check.sh %s\n' "$ROOT_ABS"
}

main "$@"
exit "$E_OK"
