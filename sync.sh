#!/usr/bin/env bash
set -e

COMMIT_TIME=$(date +"%Y-%m-%d %H:%M:%S")
BASEPATH=$(dirname "$(realpath "$0")")
cd "$BASEPATH"

# 仅在工作区有变更时才提交
if ! git diff --quiet || ! git diff --cached --quiet || [ -n "$(git ls-files --others --exclude-standard)" ]; then
    git add --all
    git commit -m "$COMMIT_TIME"
fi
if git remote | grep -q "^codeup$"; then
    git pull --rebase codeup main
    git push codeup main
fi
if git remote | grep -q "^origin$"; then
    git pull --rebase origin main
    git push origin main
fi
