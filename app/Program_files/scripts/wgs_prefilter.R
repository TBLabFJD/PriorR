# Background WGS prefilter (see app/R/wgs_prefilter.R):
#   Rscript wgs_prefilter.R <input> <profile.json> <output.tsv> <cache_dir> <progress.json> <memory_limit> <threads>
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "R", "wgs_prefilter.R"))
cfg <- jsonlite::fromJSON(args[2], simplifyVector = FALSE)
ok <- run_wgs_prefilter(args[1], cfg, args[3], args[4], args[5], args[6], as.integer(args[7]))
quit(status = if (isTRUE(ok)) 0 else 1)
