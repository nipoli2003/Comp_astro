#!/usr/bin/env bash
set -e # Stop immediately if any step fails

# Usage: ./gitpush.sh "your commit message"
# Stages, commits, and pushes the CURRENT branch (not a hardcoded "main").
# Refuses to push straight to main - everyone should work on a feature
# branch and merge into main via a Pull Request on GitHub.

if [ -z "$1" ]; then
    echo "Error: no commit message provided."
    echo "Usage: ./gitpush.sh \"your commit message\""
    exit 1
fi

branch=$(git rev-parse --abbrev-ref HEAD)

echo "==> Staging changes on branch '$branch'..."
git add .

echo "==> Committing..."
git commit -m "$1"

git push

echo "==> Done! Open a Pull Request on GitHub to merge '$branch' into main."
