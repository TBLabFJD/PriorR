# Reading of the input files. All of them return plain data frames / vectors.

# SNV file of the FJD pipeline (also the output of the WGS prefilter), with the column names of
# the file (samples as "P001_GT")
read_variants <- function(path) {
  data.table::fread(path, quote = "", header = TRUE, sep = "\t", na.strings = c("", ".", "NA"),
                    check.names = FALSE, data.table = FALSE)
}

# CNV file of AnnotSV. Split rows give the class of their CNV as "full=5". The genotype of each
# sample column (GT:FT:...) is added as <sample>_GT, like in the SNV files.
read_cnvs <- function(path) {
  df <- data.table::fread(path, quote = "", header = TRUE, sep = "\t", na.strings = c("", ".", "NA"),
                          check.names = FALSE, data.table = FALSE)
  if ("ACMG_class" %in% names(df)) df$ACMG_class <- gsub("full=", "", as.character(df$ACMG_class))
  for (s in cnv_samples(df)) {
    if (!paste0(s, "_GT") %in% names(df)) df[[paste0(s, "_GT")]] <- sub(":.*$", "", as.character(df[[s]]))
  }
  df
}

# Sample columns of an AnnotSV file: the ones named as the samples of Samples_ID
cnv_samples <- function(df) {
  if (!"Samples_ID" %in% names(df)) return(character(0))
  ids <- unique(unlist(strsplit(stats::na.omit(as.character(df$Samples_ID)), ",", fixed = TRUE)))
  intersect(trimws(ids), names(df))
}

# Gene list: one gene per line (or first column of a csv/tsv)
read_gene_list <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  genes <- trimws(sub("[,\t;].*$", "", lines))
  unique(genes[genes != "" & !startsWith(genes, "#")])
}

# BED regions: chrom, start (0-based), end
read_bed <- function(path) {
  bed <- data.table::fread(path, header = FALSE, sep = "\t", select = 1:3, colClasses = list(character = 1),
                           data.table = FALSE)
  bed <- bed[!grepl("^(track|browser|#)", bed[[1]]), ]
  names(bed) <- c("chrom", "start", "end")
  bed
}

# PED file (family, id, father, mother, sex, phenotype), tab or space separated
read_ped <- function(path) {
  utils::read.table(path, header = FALSE, colClasses = "character", comment.char = "#")
}

# Variant ids of a VCF, plain or .gz (multiallelic sites split), comparable with variant_ids().
# DuckDB reads the compressed file itself (no zcat: works on Windows too)
read_vcf_ids <- function(path, name = path) {
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  vcf <- DBI::dbGetQuery(con, sprintf(
    "SELECT column0 AS CHROM, column1 AS POS, column3 AS REF, column4 AS ALT
     FROM read_csv('%s', delim = '\t', header = false, comment = '#', all_varchar = true, null_padding = true,
                   compression = '%s')", gsub("'", "''", path), if (grepl("\\.gz$", name)) "gzip" else "none"))
  alts <- strsplit(vcf$ALT, ",", fixed = TRUE)
  n <- lengths(alts)
  variant_ids(rep(vcf$CHROM, n), rep(vcf$POS, n), rep(vcf$REF, n), unlist(alts))
}

# Small gzipped TSV tables shipped with PriorR (base R: no extra package)
read_tsv_gz <- function(path) {
  utils::read.delim(gzfile(path), sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}
