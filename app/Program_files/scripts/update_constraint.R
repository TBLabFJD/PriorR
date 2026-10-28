# Gene constraint of gnomAD v4.1 (LOEUF, pLI, missense Z) for PriorR: one row per gene (MANE Select
# transcript, or the canonical one). Run from the app/ folder (needs internet, ~100 MB):
#   Rscript Program_files/scripts/update_constraint.R [gnomad.v4.1.constraint_metrics.tsv]
args <- commandArgs(trailingOnly = TRUE)
url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
src <- if (length(args)) args[1] else {
  f <- tempfile(fileext = ".tsv")
  utils::download.file(url, f, mode = "wb", quiet = TRUE)
  f
}
x <- data.table::fread(src, select = c("gene", "transcript", "canonical", "mane_select", "lof.pLI", "lof.oe_ci.upper", "mis.z_score"),
                       data.table = FALSE)
x <- x[grepl("^ENST", x$transcript) & !is.na(x$gene) & x$gene != "", ]
# MANE Select first, then canonical
x <- x[order(x$gene, -(x$mane_select %in% c(TRUE, "true")), -(x$canonical %in% c(TRUE, "true"))), ]
x <- x[!duplicated(x$gene), ]
out <- data.frame(SYMBOL = x$gene, gnomAD_LOEUF = round(x$lof.oe_ci.upper, 3), gnomAD_pLI = signif(x$lof.pLI, 3),
                  gnomAD_mis_z = round(x$mis.z_score, 2))
dir <- file.path("Program_files", "Constraint")
dir.create(dir, showWarnings = FALSE)
con <- gzfile(file.path(dir, "gnomad_v4.1_constraint.tsv.gz"), "w")
utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(con)
message(nrow(out), " genes")
