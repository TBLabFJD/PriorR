# Downloads the exon coordinates of the transcripts (RefSeq and GENCODE tables of the UCSC Genome
# Browser, ~30 MB) into the shared folder: shared/annotations/transcripts_<build>.parquet.
# PriorR runs it in the background the first time a gene structure is drawn. By hand (e.g. Docker):
#   Rscript Program_files/scripts/get_transcripts.R [GRCh38 | GRCh37]
args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
setwd(file.path(script_dir, "..", ".."))   # the app folder, as when PriorR runs
for (f in c("config.R", "links.R", "wgs_prefilter.R", "projects.R", "gene_structure.R")) source(file.path("R", f))
build <- if (length(args)) args[1] else genome_build()
db <- transcript_db_path(build)
lock <- paste0(db, ".downloading")
dir.create(dirname(db), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(list(pid = Sys.getpid(), time = format(Sys.time())), auto_unbox = TRUE), lock)
ok <- tryCatch({
  build_transcript_db(build, db)
  unlink(paste0(db, ".error"))
  cat("Transcripts of", build, "in", db, "\n")
  TRUE
}, error = function(e) {
  writeLines(conditionMessage(e), paste0(db, ".error"))
  message("Error: ", conditionMessage(e))
  FALSE
})
unlink(lock)
quit(status = if (ok) 0 else 1)
