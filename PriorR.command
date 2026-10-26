#!/bin/bash
# Starts PriorR (macOS; double click). Needs R installed: https://cran.r-project.org/bin/macosx/
cd "$(dirname "$0")" || exit 1
RSCRIPT=$(command -v Rscript || echo /usr/local/bin/Rscript)
if [ ! -x "$RSCRIPT" ]; then
  echo "R is not installed. Install it from https://cran.r-project.org/bin/macosx/"
  read -r -p "Press Enter to close"
  exit 1
fi
exec "$RSCRIPT" launch.R
