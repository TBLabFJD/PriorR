# Runs GLOWgenes with Docker (image tblabfjd/glowgenes) in the background for PriorR:
#   Rscript run_glowgenes.R <workdir with genes.txt> <ranking file to write> <image>
# Writes <workdir>/status.json (stage: running / done / error) and the ranking (SYMBOL, rank: panel
# genes 0, then the GLOWgenes prioritization).
args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
workdir <- normalizePath(args[1]); out_file <- args[2]; image <- args[3]
# keeps who started the run (written by PriorR)
started <- tryCatch(jsonlite::fromJSON(file.path(workdir, "status.json")), error = function(e) list())
status <- function(stage, message = "") jsonlite::write_json(list(stage = stage, message = message, user = started$user %||% "",
                                                                  time = started$time %||% format(Sys.time())),
                                                            file.path(workdir, "status.json"), auto_unbox = TRUE)
status("running", "Running GLOWgenes (it can take a long time)...")
script <- paste(
  "NETS=$(dirname $(find /opt/GLOWgenesNets -name pathwayREACTOME_HGNCnets.txt | head -1))",
  "&& sed \"s#/PATH#$NETS#\" /opt/GLOWgenes/networks_knowledgeCategories.cfg > /data/networks.cfg",
  "&& mkdir -p /data/out && python /opt/GLOWgenes/GLOWgenes.py -i /data/genes.txt -n /data/networks.cfg -o /data/out")
# Linux/macOS: the files are written with the user's id, not root
user <- if (.Platform$OS.type == "unix") c("--user", paste0(system("id -u", intern = TRUE), ":", system("id -g", intern = TRUE)), "-e", "HOME=/data")
code <- system2("docker", c("run", "--rm", user, "-v", shQuote(paste0(workdir, ":/data")), image, "bash", "-c", shQuote(script)),
                stdout = file.path(workdir, "glowgenes.log"), stderr = file.path(workdir, "glowgenes.log"))
result <- file.path(workdir, "out", "GLOWgenes_prioritization_Random.txt")
if (!identical(code, 0L) || !file.exists(result)) {
  log <- if (file.exists(file.path(workdir, "glowgenes.log"))) tail(readLines(file.path(workdir, "glowgenes.log")), 3) else ""
  status("error", paste(c("GLOWgenes failed:", log), collapse = " "))
  quit(status = 1)
}
scripts_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
data_dir <- dirname(scripts_dir)   # Program_files, used by glowgenes.R
source(file.path(scripts_dir, "..", "..", "R", "glowgenes.R"))
genes <- readLines(file.path(workdir, "genes.txt"))
ranking <- glow_output_ranking(result, genes)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(ranking, paste0(out_file, ".tmp"), sep = "\t")
file.rename(paste0(out_file, ".tmp"), out_file)
status("done", sprintf("GLOWgenes ranking with %d genes", nrow(ranking)))
# PriorR's own work folder (glowgenes/work/<name>) is removed once the ranking is saved
if (basename(dirname(workdir)) == "work") unlink(workdir, recursive = TRUE)
