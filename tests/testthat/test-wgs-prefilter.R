run_prefilter <- function(df, profile) {
  dir <- tempfile("wgs-"); dir.create(dir)
  input <- file.path(dir, "wgs.tsv"); data.table::fwrite(df, input, sep = "\t", na = "")
  ok <- run_wgs_prefilter(input, profile, file.path(dir, "out.tsv"), file.path(dir, "cache"), file.path(dir, "progress.json"))
  list(ok = ok, out = if (file.exists(file.path(dir, "out.tsv"))) read_variants(file.path(dir, "out.tsv")),
       progress = jsonlite::fromJSON(file.path(dir, "progress.json")))
}

test_that("the WGS prefilter keeps candidates and rescues ClinVar / SpliceAI variants", {
  df <- make_variants(6, Genomic_region = c("EXONIC", "INTRONIC", "INTRONIC", "INTRONIC", "EXONIC", "SPLICING"),
                      gnomADg_AF_popmax = c(0.001, 0.001, 0.001, 0.5, 0.5, NA),
                      CLNSIG = c(NA, "Pathogenic", NA, NA, NA, "Benign"))
  df$SpliceAI_DS_Max <- c(NA, NA, 0.8, 0.9, NA, NA)
  profile <- list(regions = list("EXONIC", "SPLICING"), canonical_only = FALSE,
                  max_af = list(gnomADg_AF_popmax = 0.01), min_cadd = NULL, exclude_clinvar_benign = TRUE,
                  rescue_clinvar_pathogenic = TRUE, rescue_spliceai_min = 0.5, genes = list(), max_rows = 1000)
  r <- run_prefilter(df, profile)
  expect_true(r$ok)
  # 1 exonic rare; 2 ClinVar pathogenic intronic; 3 SpliceAI rare; not 4 (common), 5 (common), 6 (benign)
  expect_equal(r$out$POS, df$POS[1:3])
  expect_equal(r$progress$stage, "done")
  expect_equal(r$progress$total_rows, 6)
  # same columns as the input, so PriorR reads it as a normal SNV file
  expect_equal(names(r$out), names(df))
})

test_that("more rows than the limit: nothing is loaded and the error says why (never a partial table)", {
  df <- make_variants(5, CLNSIG = c(NA, NA, "Pathogenic", NA, NA), CADD_PHRED = c(10, 40, 5, 30, 20))
  profile <- list(regions = list("EXONIC"), max_af = list(), rescue_clinvar_pathogenic = TRUE, max_rows = 2)
  r <- run_prefilter(df, profile)
  expect_equal(r$progress$stage, "error")
  expect_match(r$progress$message, "5 rows pass the load filters, more than the safety limit of 2 rows")
  expect_null(r$out)
  # within the limit: everything
  r <- run_prefilter(df, modifyList(profile, list(max_rows = 5)))
  expect_equal(nrow(r$out), 5)
})

test_that("errors are reported in the progress file", {
  dir <- tempfile("wgs-"); dir.create(dir)
  ok <- run_wgs_prefilter(file.path(dir, "missing.tsv"), list(), file.path(dir, "out.tsv"), dir, file.path(dir, "p.json"))
  expect_false(ok)
  expect_equal(jsonlite::fromJSON(file.path(dir, "p.json"))$stage, "error")
})

test_that("the load keeps the genes of the panels and drops the excluded ones", {
  df <- make_variants(4, SYMBOL = c("A", "B", "C", NA))
  r <- run_prefilter(df, list(genes = list("A", "B", "C"), exclude_genes = list("B")))
  expect_equal(r$out$SYMBOL, c("A", "C"))
  r <- run_prefilter(df, list(exclude_genes = list("B")))
  expect_equal(r$out$POS, df$POS[c(1, 3, 4)])        # variants without gene are kept
})

test_that("a file that is not a variant table gives a clear error", {
  dir <- tempfile("wgs-"); dir.create(dir)
  f <- file.path(dir, "notes.txt"); writeLines(c("just some text", "not a table"), f)
  ok <- run_wgs_prefilter(f, list(), file.path(dir, "out.tsv"), dir, file.path(dir, "p.json"))
  expect_false(ok)
  expect_match(jsonlite::fromJSON(file.path(dir, "p.json"))$message, "can't be read as a variant table")
})
