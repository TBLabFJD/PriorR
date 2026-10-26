#!/bin/bash
# Starts PriorR (Linux). Needs R installed (e.g. sudo apt install r-base).
cd "$(dirname "$0")" || exit 1
if ! command -v Rscript >/dev/null; then
  echo "R is not installed (e.g. sudo apt install r-base)"
  exit 1
fi
exec Rscript launch.R
