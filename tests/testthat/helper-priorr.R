# Loads the app code with a temporary database and builds synthetic test data
# (no patient data in the repository).

library(shiny)

app_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "app"))
test_files_dir <- tempfile("priorr-files-")
dir.create(test_files_dir)
Sys.setenv(PRIORR_DIR = file.path(app_dir, "Program_files"),
           PRIORR_SHARED_DIR = tempfile("priorr-shared-"),
           PRIORR_CACHE_DIR = tempfile("priorr-cache-"),
           PRIORR_DATA_DIRS = paste0("Tests=", test_files_dir),
           # nothing is written in the user's folders while testing
           PRIORR_USER_CONFIG = tempfile("priorr-cfg-"), PRIORR_LOGIN_DIR = tempfile("priorr-login-"),
           PRIORR_LOCAL_SHARED = tempfile("priorr-local-"))
for (f in list.files(file.path(app_dir, "R"), full.names = TRUE)) source(f, local = FALSE)

# A variant table like the ones of the FJD pipeline. Each argument overrides a column.
make_variants <- function(n = 6, ...) {
  df <- data.frame(
    CHROM = "chr1", POS = seq(1000, by = 100, length.out = n), REF = "A", ALT = "G",
    SYMBOL = "GENE1", VARIANT_CLASS = "SNV", Genomic_region = "EXONIC", CANONICAL = "YES",
    Feature = "NM_000001.1", Consequence = "missense_variant", CLNSIG = NA_character_,
    OMIM_phenotype = "Disease, Autosomal recessive", gnomADg_AF = 0.0001, gnomADg_AF_popmax = 0.0002,
    gnomADe_AF = 0.0001, FJD_MAF_AF = 0.0001, CADD_PHRED = 20, MutScore = 0.5,
    N_Pathogenic_pred = 5, N_Benign_pred = 5, N_predictions = 10, ada_score = NA_real_, rf_score = NA_real_,
    GLOWgenes = 10, `P001_GT` = "0/1", `P001_DP` = 30,
    stringsAsFactors = FALSE, check.names = FALSE)
  overrides <- list(...)
  for (col in names(overrides)) df[[col]] <- rep_len(overrides[[col]], n)
  df$variant_id <- paste(df$CHROM, df$POS, df$REF, df$ALT, sep = "_")
  df
}

# A fresh shared folder for each test that needs one
fresh_store <- function() {
  shared_dir <<- tempfile("priorr-shared-")
  dir.create(shared_dir)
  invisible(shared_dir)
}

# Value of an "Open file" input choosing a file of the test folder
open_file_value <- function(df, ext = ".tsv") {
  f <- tempfile(tmpdir = test_files_dir, fileext = ext)
  data.table::fwrite(df, f, sep = "\t", na = "")
  list(root = "Tests", files = list(`0` = list("", basename(f))))
}
