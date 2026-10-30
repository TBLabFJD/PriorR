# Projects: other samples of the project with each variant (projects.R), synthetic files

write_sample <- function(name, df) {
  path <- file.path(test_files_dir, name)
  data.table::fwrite(df, path, sep = "\t", na = "")
  path
}

test_that("a project indexes its files and counts each variant in the other samples", {
  fresh_store()
  # three files: S1 (the sample analysed), S2 and a trio file with S3, S4, S5 (S5 without the variants)
  s1 <- write_sample("S1.tsv", make_variants(3, POS = c(1000, 2000, 100000000), P001_GT = c("0/1", "1/1", "0/1")) |>
                       stats::setNames(sub("P001", "S1", names(make_variants(1)))))
  s2 <- write_sample("S2.tsv", make_variants(2, POS = c(1000, 2000), P001_GT = c("1/1", "0/1")) |>
                       stats::setNames(sub("P001", "S2", names(make_variants(1)))))
  trio <- make_variants(3, POS = c(1000, 1000, 100000000))
  trio$Feature <- c("NM_1", "NM_2", "NM_1")   # two transcripts of the same variant count once
  names(trio) <- sub("P001", "S3", names(trio))
  trio$S4_GT <- c("0/1", "0/1", "0|1"); trio$S5_GT <- c("0/0", "0/0", "./.")
  trio$S3_DV_GT <- "0/1"   # a caller of S3, not a sample
  s3 <- write_sample("trio.tsv", trio)

  id <- create_project("Retina 2026", "ana")
  expect_error(create_project("retina 2026", "ana"), "already")
  expect_equal(unname(project_choices()), id)
  expect_equal(add_project_files(id, c(s1, s2, s3, s2), "ana"), 3)
  expect_equal(project_status(id)$indexed, c(FALSE, FALSE, FALSE))

  index_project(id, cache_dir = tempfile())
  st <- project_status(id)
  expect_true(all(st$indexed))
  expect_equal(st$samples, c("S1", "S2", "S3, S4, S5"))
  expect_equal(project_indexing(id)$stage, "done")

  # analysing S1: the other samples are S2..S5
  df <- read_variants(s1)
  counts <- project_counts(id, df, exclude = "S1")
  expect_equal(counts$total, 4)
  c1 <- counts$counts[match(variant_key("chr1", c("1000", "2000", "100000000"), "A", "G"), counts$counts$key), ]
  expect_equal(c1$n, c(3L, 1L, 2L))      # 1000: S2, S3, S4 · 2000: S2 · 1e8: S3, S4 (POS without scientific notation)
  expect_equal(c1$hom, c(1L, 0L, 0L))

  out <- filter_snv(df, snv_params(project = counts))
  expect_equal(out$Project_samples, c(3L, 1L, 2L))
  expect_equal(out$Project_hom, c(1L, 0L, 0L))
  expect_equal(attr(out, "display_columns")[match("SYMBOL", attr(out, "display_columns")) + 1:3],
               c("Variants_in_gene", "Project_samples", "Project_hom"))
  expect_equal(filter_snv(df, snv_params(project = counts, max_project = 1))$POS, 2000)

  who <- project_variant_samples(id, "chr1", 1000, "A", "G")
  expect_equal(who$sample, c("S1", "S2", "S3", "S4"))
  expect_equal(who$gt, c("het", "hom", "het", "het"))

  # a changed file is indexed again; a removed one leaves the counts
  Sys.sleep(1.1)
  write_sample("S2.tsv", stats::setNames(make_variants(1, POS = 5000), sub("P001", "S2", names(make_variants(1)))))
  expect_false(project_status(id)$indexed[2])
  remove_project_file(id, s3)
  index_project(id, cache_dir = tempfile())
  expect_equal(project_counts(id, df, exclude = "S1")$counts$n, integer(0))
  delete_project(id)
  expect_equal(length(list_projects()), 0)
})

test_that("genotypes: het, hom, hemizygous, no alternative allele", {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  gts <- c("0/1", "1/1", "1|1", "0|2", "1", "0/0", "./.", "./1", ".")
  res <- DBI::dbGetQuery(con, sprintf("SELECT %s AS gt FROM (SELECT unnest([%s]) AS x_GT)", gt_sql("x_GT"),
                                      paste(sql_str(gts), collapse = ", ")))$gt
  expect_equal(res, c("het", "hom", "hom", "het", "hom", NA, NA, "het", NA))
})
