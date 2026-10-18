#!/usr/bin/env bash
# ============================================================================
# repro_check.sh —— 文献复现项目一致性自检
# 上位规范：REPRO_SPEC.md（通用规范 v2.1，L0）§11.4
#
# 用法：
#   tools/repro_check.sh <项目目录>        # 例：tools/repro_check.sh 20261016_GBM_lactate_NFATc4_SPP1_macrophage
#   tools/repro_check.sh .                 # 在项目目录内执行
#   tools/repro_check.sh --all             # 扫描仓库下所有项目目录
#
# 退出码：0 = 全部通过；1 = 存在需修正的问题
# ============================================================================
set -uo pipefail

FAIL=0
note() { printf '  %s\n' "$*"; }
sec()  { printf '\n[%s] %s\n' "$1" "$2"; }
bad()  { printf '  ❌ %s\n' "$*"; FAIL=1; }
warn() { printf '  ⚠️  %s\n' "$*"; }
ok()   { printf '  ✅ %s\n' "$*"; }

check_project() {
  local dir="$1"
  printf '\n════════════════════════════════════════════════════════════\n'
  printf '项目：%s\n' "$dir"
  printf '════════════════════════════════════════════════════════════\n'
  [ -d "$dir" ] || { bad "目录不存在：$dir"; return; }
  cd "$dir" || return

  # --- 1. 必需结构 -----------------------------------------------------------
  sec 1 "必需文件与目录"
  [ -f PROJECT_SPEC.md ] && ok "PROJECT_SPEC.md 存在（L1 项目规范）" || bad "缺少 PROJECT_SPEC.md（见 L0 附录 A 模板）"
  ls 00_*.md >/dev/null 2>&1 && ok "00 号总览存在" || bad "缺少 00_*.md 总览文档"
  for d in raw processed results figures code logs; do
    [ -d "$d" ] || warn "目录未创建：$d/（若尚未执行代码可忽略）"
  done
  [ -f logs/DECISIONS.md ] && ok "logs/DECISIONS.md 存在" || warn "缺少 logs/DECISIONS.md（L0 §10.2 要求运行时创建）"

  # --- 2. 代码围栏配对 -------------------------------------------------------
  sec 2 "Markdown 代码围栏配对（必须为偶数）"
  local odd=0
  for f in *.md; do
    [ -e "$f" ] || continue
    local n; n=$(grep -c '^```' "$f" 2>/dev/null | head -1 || true); n=${n:-0}
    if [ $((n % 2)) -ne 0 ]; then bad "$f：围栏数 $n（奇数，未闭合）"; odd=1; fi
  done
  [ "$odd" -eq 0 ] && ok "全部配对"

  # --- 3. 交叉引用有效性 -----------------------------------------------------
  sec 3 "文档交叉引用（NN_*.md 必须指向真实文件）"
  local miss=0 ref
  while read -r ref; do
    [ -z "$ref" ] && continue
    case "$ref" in */*|*"..."*|*"{"*) continue;; esac          # 跳过路径/省略号/占位写法
    if [ ! -f "$ref" ]; then bad "断链：$ref"; miss=1; fi
  done < <(grep -ohP '(?<![A-Za-z0-9_/])[0-9][0-9]_[^ )）`"]*\.md' *.md 2>/dev/null | sort -u)
  [ "$miss" -eq 0 ] && ok "无断链"

  # --- 4. 占位符 / 未完成标记 ------------------------------------------------
  sec 4 "占位符与未完成标记（应为 0）"
  local ph
  ph=$(grep -nE '(^|[^A-Za-z])TODO([^A-Za-z]|$)|待补|待填|此处省略|FIXME' *.md 2>/dev/null \
        | grep -v '无 TODO\|禁止 TODO' || true)
  if [ -n "$ph" ]; then bad "发现未完成标记："; printf '%s\n' "$ph" | head -5; else ok "无"; fi

  # --- 5. 文件名非法字符（仅项目产物目录与顶层）------------------------------
  sec 5 "文件名非法字符（< > { } : * ? | 空格）"
  local badname
  badname=$(ls -1 | grep -E '[<>{}:*?|]' || true)               # 顶层：文档允许中文，但不应含这些符号
  local b2; b2=$(find raw processed results figures code logs -mindepth 1 -maxdepth 2 2>/dev/null \
                 | grep -E '[<>{}:*?| ]' || true)
  if [ -n "$badname$b2" ]; then bad "非法文件名："; printf '%s\n%s\n' "$badname" "$b2" | sed '/^$/d' | head -8
  else ok "无（顶层文档的空白/中文按 L0 §3.2 豁免）"; fi

  # --- 6. 头部四要素 ---------------------------------------------------------
  sec 6 "模块文档头部四要素（目标文献/模块对应原文/前置依赖/复现难度）"
  local inc=0
  for f in [0-9][0-9]_*.md; do
    [ -e "$f" ] || continue
    case "$f" in 00_*) continue;; esac      # 00 为总览/索引，不适用模块文档四要素
    local c; c=$(grep -c '目标文献\|模块对应原文\|前置依赖\|复现难度' "$f")
    if [ "$c" -lt 4 ]; then bad "$f：仅命中 $c/4"; inc=1; fi
  done
  [ "$inc" -eq 0 ] && ok "全部齐全"

  # --- 7. 已知旧名残留（L0/L1 通用反模式）------------------------------------
  sec 7 "旧路径/旧名残留（data/raw、data/processed）"
  local legacy
  legacy=$(grep -rn 'data/raw/\|data/processed/' --exclude=PROJECT_SPEC.md *.md 2>/dev/null \
           | grep -v '不使用\|指同一处\|均指\|早期草稿\|已废弃\|拒绝\|旧写法' || true)
  if [ -n "$legacy" ]; then warn "存在历史写法（确认是否为说明性引用）："; printf '%s\n' "$legacy" | head -5
  else ok "无"; fi

  # --- 8. 缺失项显式标注 -----------------------------------------------------
  sec 8 "文献缺失项是否显式标注（应为正数）"
  local cnt; cnt=$(grep -oh '原文未报告\|文献未报告\|原文未写' *.md 2>/dev/null | wc -l)
  if [ "$cnt" -gt 0 ]; then ok "$cnt 处（说明缺失项已被显式处理）"; else warn "0 处 —— 请确认是否真的没有缺失项"; fi

  cd - >/dev/null || true
}

# ---------------------------------------------------------------------------
if [ $# -eq 0 ]; then
  echo "用法: $0 <项目目录> | --all"; exit 2
fi

if [ "$1" = "--all" ]; then
  for d in */; do
    d="${d%/}"
    [ -f "$d/PROJECT_SPEC.md" ] && check_project "$d"
  done
else
  check_project "$1"
fi

printf '\n════════════════════════════════════════════════════════════\n'
if [ "$FAIL" -eq 0 ]; then
  printf '结果：✅ 通过（无阻断性问题）\n'; exit 0
else
  printf '结果：❌ 存在需修正的问题（见上方 ❌ 标记）\n'; exit 1
fi
