# Background indexing of the files of a project (see app/R/projects.R):
#   Rscript index_project.R <project id>
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
setwd(file.path(script_dir, "..", ".."))   # the app folder, as when PriorR runs
for (f in c("config.R", "store.R", "wgs_prefilter.R", "inheritance.R", "projects.R")) source(file.path("R", f))
ok <- index_project(args[1], wgs_cache_dir, wgs_memory_limit, as.integer(wgs_threads))
quit(status = if (isTRUE(ok)) 0 else 1)
