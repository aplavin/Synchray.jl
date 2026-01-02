#!/usr/bin/env bash
set -e

tmp="$(mktemp -d)"
git clone --no-local --no-hardlinks "$(pwd)" "$tmp/Synchray.jl" >/dev/null
cd "$tmp/Synchray.jl"

# Rewrite history:
# - Remove the specified files from all commits
# - Round (floor) author/committer dates to the first day of the month at 00:00:00 UTC
#
# This rewrites commit hashes; pushing requires --force.

git filter-branch --force --prune-empty --tag-name-filter cat \
  --index-filter 'git rm -r --cached --ignore-unmatch DESIGN.md PLAN.md USECASES.md publish.sh' \
    --env-filter '
ts=$(git show -s --format=%at "$GIT_COMMIT")
ym=$(date -u -r "$ts" +%Y-%m)
export GIT_AUTHOR_DATE=${ym}-01T00:00:00Z

ts=$(git show -s --format=%ct "$GIT_COMMIT")
ym=$(date -u -r "$ts" +%Y-%m)
export GIT_COMMITTER_DATE=${ym}-01T00:00:00Z
' \
  -- --all >/dev/null

# Cleanup filter-branch backup refs and aggressively gc to shrink the repo.
rm -rf .git/refs/original .git/logs/refs/original
git reflog expire --expire=now --all >/dev/null
git gc --prune=now --aggressive >/dev/null

git remote set-url origin https://github.com/JuliaAPlavin/Synchray.jl || git remote add origin https://github.com/JuliaAPlavin/Synchray.jl
git push --force --all origin
git push --force --tags origin
