test_that("ClinVar pathogenic variants are red", {
  df <- make_variants(3, CLNSIG = c("Pathogenic", "Likely_pathogenic", "Benign"))
  expect_equal(priority_class(df), c("red", "red", NA))
})

test_that("rare, predicted pathogenic, high CADD and MutScore is red", {
  df <- make_variants(2, N_Pathogenic_pred = 8, CADD_PHRED = c(30, 20), MutScore = 0.9, gnomADg_AF_popmax = 0.0001)
  expect_equal(priority_class(df), c("red", NA))
})

test_that("low frequency (0.001-0.01) candidates are orange", {
  df <- make_variants(1, gnomADg_AF_popmax = 0.005, N_Pathogenic_pred = 8, CADD_PHRED = 30, MutScore = 0.9)
  expect_equal(priority_class(df), "orange")
})

test_that("both ClinVar names of conflicting classifications are not orange", {
  df <- make_variants(3, CLNSIG = c("Conflicting_interpretations_of_pathogenicity",
                                    "Conflicting_classifications_of_pathogenicity",
                                    "Pathogenic|risk_factor"))
  expect_equal(priority_class(df), c(NA, NA, "orange"))
})

test_that("a row matching red and orange is only red and never duplicated", {
  df <- make_variants(3, CLNSIG = c(NA, "Pathogenic", NA), gnomADg_AF_popmax = c(0.5, 0.005, 0.5),
                      N_Pathogenic_pred = 8, CADD_PHRED = 30, MutScore = 0.9)
  expect_equal(priority_class(df), c(NA, "red", NA))
  out <- prioritize(df)
  expect_equal(nrow(out), 3)
  expect_equal(out$POS[1], df$POS[2])
})

test_that("prioritization keeps the order inside each class", {
  df <- make_variants(4, CLNSIG = c(NA, "Pathogenic", NA, "Pathogenic"))
  expect_equal(prioritize(df)$POS, df$POS[c(2, 4, 1, 3)])
})

test_that("Orphanet panels rank the panel genes first", {
  panels <- orphanet_panels()
  expect_true("retinal_dystrophy_panel_glow" %in% panels)
  rank <- orphanet_rank("retinal_dystrophy_panel_glow")
  expect_equal(rank$Orphanet_GLOW_rank[rank$SYMBOL == "ABCA4"], 0)
  expect_gt(rank$Orphanet_GLOW_rank[rank$SYMBOL == "TNF"], 0)
})

test_that("GLOWgenes view: rare variants with evidence, one row per variant, ordered by the rank", {
  df <- make_variants(7, SYMBOL = c("A", "A", "B", "C", "D", "E", "F"), POS = c(1, 1, 2, 3, 4, 5, 6),
                      CANONICAL = c("", "YES", "YES", "YES", "YES", "YES", "YES"),
                      GLOWgenes = c(5, 5, 0, 2, 1, NA, 3), CADD_PHRED = c(30, 30, 25, 10, 22, 28, 30),
                      gnomADg_AF_popmax = c(0.001, 0.001, 0.0001, 0.0001, 0.0001, 0.0001, 0.2),
                      CLNSIG = c(NA, NA, NA, NA, "Benign", NA, NA))
  out <- glow_view(df)
  # C: no evidence (CADD 10); D: benign; F: common; one row of A (canonical)
  expect_equal(out$SYMBOL, c("B", "A", "E"))
  expect_equal(out$GLOW_rank, c(0, 5, NA))
  expect_equal(out$CANONICAL[out$SYMBOL == "A"], "YES")
  expect_true(has_glow_column(df))
  # ranking from an Orphanet panel: its genes first
  df2 <- make_variants(2, SYMBOL = c("TNF", "ABCA4"), CADD_PHRED = 30)
  df2$GLOWgenes <- NULL
  expect_false(has_glow_column(df2))
  expect_equal(glow_view(df2, "orphanet:retinal_dystrophy_panel_glow")$SYMBOL, c("ABCA4", "TNF"))
})

test_that("GLOWgenes view through the SNV filter keeps the colors and shows the rank first", {
  df <- make_variants(3, SYMBOL = c("A", "B", "C"), GLOWgenes = c(9, 0, 4), CADD_PHRED = 30)
  out <- filter_snv(df, snv_params(columns = c("SYMBOL"), glow = list(sources = NULL)))
  expect_equal(out$SYMBOL, c("B", "C", "A"))
  expect_equal(attr(out, "display_columns")[1], "GLOW_rank")
  expect_length(attr(out, "priority"), 3)
})

test_that("GLOWgenes rankings of PanelApp (Sci Rep 2026), specific genes and SGDS", {
  choices <- glow_rank_choices()
  expect_true(all(c("Orphanet", "PanelApp (GLOWgenes, Sci Rep 2026)") %in% names(choices)))
  expect_true("panelapp:Retinal_disorders_GA" %in% unlist(choices))
  ranks <- glow_ranking("panelapp:Retinal_disorders_GA")
  expect_equal(ranks$rank[ranks$SYMBOL == "ERCC6"], 0)
  expect_equal(ranks$SYMBOL[which(ranks$rank == 1)], "FSCN2")
  df <- make_variants(2, SYMBOL = c("FSCN2", "ERCC6"))
  expect_equal(glow_rank(df, "panelapp:Retinal_disorders_GA"), c(1, 0))

  catalog <- panel_catalog()
  specific <- catalog[catalog$source == "GLOWgenes specific genes (Sci Rep 2026)", ]
  expect_gt(nrow(specific), 10)
  expect_true(any(grepl("Ciliopathies - class-specific genes", specific$name)))
  expect_equal(unname(lengths(specific$genes[specific$id == "glowspecific:class:Ciliopathies"])), 130)

  sg <- add_sgds(make_variants(2, SYMBOL = c("SNHG22", "NOT_A_GENE")))
  expect_equal(sg$SGDS, c(1, NA))
})

test_that("a GLOWgenes run becomes a ranking: panel genes 0, then its prioritization", {
  f <- tempfile()
  writeLines(c("GENEX\t0.9\t1", "GENEA\t0.8\t2", "GENEY\t0.5\t3"), f)
  r <- glow_output_ranking(f, c("GENEA", "GENEB"))
  expect_equal(r$SYMBOL, c("GENEA", "GENEB", "GENEX", "GENEY"))
  expect_equal(r$rank, c(0, 0, 1, 3))
})

test_that("the max GLOWgenes rank uses the file column or the chosen rankings, in both views", {
  df <- make_variants(3, SYMBOL = c("A", "B", "C"), GLOWgenes = c(9, 0, 4), CADD_PHRED = 30)
  expect_equal(filter_snv(df, snv_params(glow_max = 4))$SYMBOL, c("B", "C"))
  expect_equal(filter_snv(df, snv_params(glow_max = 4, glow = list(sources = NULL)))$SYMBOL, c("B", "C"))
  df2 <- make_variants(2, SYMBOL = c("FSCN2", "TNF"))
  expect_equal(filter_snv(df2, snv_params(glow_max = 1, glow_sources = "panelapp:Retinal_disorders_GA"))$SYMBOL, "FSCN2")
})

test_that("GLOWgenes runs in progress or failed are listed from the shared folder", {
  for (r in list(list("Retina", "running"), list("Old", "error"), list("Ok", "done"))) {
    dir.create(file.path(glow_runs_dir(), "work", r[[1]]), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(list(stage = r[[2]], message = "m", user = "pablo", time = "2026-10-23 16:00:00"),
                         file.path(glow_runs_dir(), "work", r[[1]], "status.json"), auto_unbox = TRUE)
  }
  runs <- glow_run_statuses()
  expect_setequal(runs$name, c("Retina", "Old"))
  expect_equal(runs$user[runs$name == "Retina"], "pablo")
  unlink(file.path(glow_runs_dir(), "work"), recursive = TRUE)
  expect_equal(nrow(glow_run_statuses()), 0)
})

test_that("GLOWgenes view can hide the genes of the panel (rank 0)", {
  df <- make_variants(3, SYMBOL = c("A", "B", "C"), GLOWgenes = c(9, 0, 4), CADD_PHRED = 30)
  out <- filter_snv(df, snv_params(glow = list(sources = NULL, hide_rank0 = TRUE)))
  expect_equal(out$SYMBOL, c("C", "A"))
})
