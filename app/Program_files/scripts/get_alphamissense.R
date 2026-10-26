# AlphaMissense predictions for PriorR (Cheng et al., Science 2023; CC BY 4.0, https://zenodo.org/records/10813168).
# Downloads AlphaMissense_hg38.tsv.gz (~640 MB) and writes a Parquet file (CHROM, POS, REF, ALT, score, class)
# to the shared folder (setting alphamissense_file), used when a variant file has no AlphaMissense column.
# Run from the app/ folder (needs internet the first time):
#   Rscript Program_files/scripts/get_alphamissense.R [AlphaMissense_hg38.tsv.gz]
suppressMessages(library(shiny))
for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else {
  f <- tempfile(fileext = ".tsv.gz")
  options(timeout = 3600)
  utils::download.file("https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz?download=1", f, mode = "wb", quiet = TRUE)
  f
}
dir.create(dirname(alphamissense_file), showWarnings = FALSE, recursive = TRUE)
tmp <- paste0(alphamissense_file, ".tmp")
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
tmp_dir <- file.path(wgs_cache_dir, "duckdb_tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
invisible(DBI::dbExecute(con, sprintf("SET memory_limit = '%s'; SET threads = %s; SET temp_directory = '%s';",
                            wgs_memory_limit, wgs_threads, gsub("'", "''", tmp_dir))))
# 3 lines of license comments, then the header "#CHROM POS REF ALT genome uniprot_id ..."
invisible(DBI::dbExecute(con, sprintf(
  "COPY (SELECT \"#CHROM\" AS CHROM, CAST(POS AS INTEGER) AS POS, REF, ALT,
          CAST(am_pathogenicity AS FLOAT) AS am_pathogenicity, am_class
        FROM read_csv('%s', delim = '\t', skip = 3, header = true, all_varchar = true)
        ORDER BY CHROM, POS)
   TO '%s' (FORMAT parquet, COMPRESSION zstd)", src, tmp)))
DBI::dbDisconnect(con, shutdown = TRUE)
invisible(file.rename(tmp, alphamissense_file))
message("AlphaMissense: ", alphamissense_file)
