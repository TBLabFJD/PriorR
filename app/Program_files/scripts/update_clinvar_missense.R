# Pathogenic / likely pathogenic missense variants of ClinVar (GRCh38, at least one star) by transcript
# and protein position, for the ACMG criteria PS1 (same amino acid change) and PM5 (other pathogenic
# change at the same residue) in acmg.R. Run from the app/ folder (needs internet, ~450 MB):
#   Rscript Program_files/scripts/update_clinvar_missense.R [variant_summary.txt.gz]
args <- commandArgs(trailingOnly = TRUE)
url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
src <- if (length(args)) args[1] else {
  f <- tempfile(fileext = ".txt.gz")
  options(timeout = 3600)
  utils::download.file(url, f, mode = "wb", quiet = TRUE)
  f
}
cols <- c("Type", "Name", "GeneSymbol", "ClinicalSignificance", "Assembly", "Chromosome", "ReviewStatus",
          "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")
# DuckDB reads the .gz itself (no zcat: works on Windows too)
con <- DBI::dbConnect(duckdb::duckdb())
x <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM read_csv('%s', delim = '\t', header = true, quote = '', all_varchar = true)",
                                  paste(sprintf('"%s"', cols), collapse = ", "), gsub("'", "''", src)))
DBI::dbDisconnect(con, shutdown = TRUE)
x <- x[x$Assembly == "GRCh38" & x$Type == "single nucleotide variant" &
         x$ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") &
         !x$ReviewStatus %in% c("no assertion criteria provided", "no classification provided", "no classifications from unflagged records"), ]
# NM_000350.3(ABCA4):c.5882G>A (p.Gly1961Glu)
m <- regmatches(x$Name, regexec("^(NM_[0-9]+)\\.[0-9]+\\([^)]*\\):c\\.[^ ]+ \\(p\\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})\\)$", x$Name))
ok <- lengths(m) == 5
x <- x[ok, ]; m <- do.call(rbind, m[ok])
out <- data.frame(transcript = m[, 2], aa_pos = as.integer(m[, 4]), ref_aa = m[, 3], alt_aa = m[, 5], gene = x$GeneSymbol,
                  chrom = paste0("chr", x$Chromosome), pos = x$PositionVCF, ref = x$ReferenceAlleleVCF, alt = x$AlternateAlleleVCF,
                  clnsig = x$ClinicalSignificance, stringsAsFactors = FALSE)
# missense only: no stop, no synonymous, no start codon (p.Met1 changes are start lost)
out <- out[out$alt_aa != "Ter" & out$ref_aa != out$alt_aa & !(out$aa_pos == 1 & out$ref_aa == "Met"), ]
out <- unique(out[order(out$transcript, out$aa_pos), ])
dir <- file.path("Program_files", "ClinVar")
dir.create(dir, showWarnings = FALSE)
con <- gzfile(file.path(dir, "clinvar_missense_plp.tsv.gz"), "w")
utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
writeLines(format(Sys.Date()), file.path(dir, "clinvar_version.txt"))
message(nrow(out), " pathogenic / likely pathogenic missense variants of ClinVar")
