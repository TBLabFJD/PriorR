#!/bin/bash
# Starts PriorR (macOS; double click). Needs R installed: https://cran.r-project.org/bin/macosx/
cd "$(dirname "$0")" || exit 1
# Rscript of the PATH, of the R installer of CRAN (/usr/local/bin) or of the R framework (Apple silicon)
RSCRIPT=$(command -v Rscript)
for candidate in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript /Library/Frameworks/R.framework/Resources/bin/Rscript; do
  [ -z "$RSCRIPT" ] && [ -x "$candidate" ] && RSCRIPT=$candidate
done
if [ ! -x "$RSCRIPT" ]; then
  echo "R is not installed. Install it from https://cran.r-project.org/bin/macosx/"
  read -r -p "Press Enter to close"
  exit 1
fi
exec "$RSCRIPT" launch.R
