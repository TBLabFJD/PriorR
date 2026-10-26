test_that("built-in profiles: None filters nothing, WES and WGS", {
  p <- builtin_profiles()
  expect_equal(names(p), c("None", "WES", "WGS"))
  expect_equal(p$None$filters$min_cadd, 0)
  expect_null(p$None$filters$region)
  expect_equal(unlist(p$WES$load$regions), c("EXONIC", "SPLICING"))
  expect_equal(p$WES$load$min_cadd, 15)
  expect_true(p$WES$load$canonical_only)
  expect_true(p$WES$load$rescue_clinvar_pathogenic)
  expect_equal(p$WES$filters$max_gnomadg_popmax, 0.01)
  expect_equal(p$WES$filters$min_cadd, 0)          # the load already filtered: rescued variants stay visible
  expect_true(p$WGS$load$rescue_clinvar_pathogenic)
  expect_equal(unlist(p$WGS$load$regions), c("EXONIC", "SPLICING"))
  expect_null(p$WGS$load$max_rows)   # the row limit is a setting of the computer, not of the profile
})

test_that("each user saves, lists and deletes their own profiles", {
  fresh_store()
  save_profile("ana", list(name = "Retina", filters = list(min_cadd = 20), columns = list("SYMBOL")))
  expect_equal(names(user_profiles("ana")), "Retina")
  expect_equal(length(user_profiles("luis")), 0)
  expect_equal(user_profiles("ana")$Retina$filters$min_cadd, 20)
  expect_equal(user_profiles("ana")$Retina$filters$max_gnomadg_af, 1)     # completed with no-filter values
  expect_error(save_profile("ana", list(name = "WGS")))                    # built-in name
  expect_true("user:Retina" %in% profile_choices("ana"))
  delete_profile("ana", "Retina")
  expect_equal(length(user_profiles("ana")), 0)
})

test_that("panel library: Orphanet, saved panels and PanelApp", {
  fresh_store()
  catalog <- panel_catalog()
  expect_true(any(catalog$source == "Orphanet"))
  retina <- catalog$id[catalog$source == "Orphanet" & grepl("retinal dystrophy", catalog$name)]
  expect_true("ABCA4" %in% panel_genes(retina))
  save_custom_panel("Cardio FJD", c("MYH7", "TTN"))
  expect_error(save_custom_panel("Cardio FJD", "X"))
  expect_setequal(panel_genes("custom:Cardio FJD"), c("MYH7", "TTN"))
  if (any(startsWith(catalog$source, "PanelApp"))) {
    pa <- panel_catalog()[startsWith(panel_catalog()$source, "PanelApp"), ]
    expect_gt(nrow(pa), 100)
    expect_true(all(lengths(pa$genes) + lengths(pa$amber) + lengths(pa$red) > 0))
  }
})

test_that("RetiGene panels: retained genes green, candidates amber, clinical subsets", {
  catalog <- panel_catalog()
  skip_if_not(any(catalog$source == "RetiGene"))
  all <- catalog[catalog$id == "retigene:all", ]
  expect_equal(unname(lengths(all$genes)), 470)
  expect_equal(unname(lengths(all$amber)), 196)
  expect_true(all(c("RPGR", "USH2A") %in% panel_genes("retigene:rp")))
  expect_true("ABCA4" %in% panel_genes("retigene:md"))
  expect_gt(length(panel_genes("retigene:all", amber = TRUE)), 470)
})

test_that("the disease-specific genes (Benitez et al. 2026) have their own panel choices", {
  all <- unlist(panel_choices())
  without <- unlist(panel_choices(specific = FALSE))
  only <- panel_choices(specific = TRUE)
  expect_false(any(startsWith(without, "glowspecific:")))
  expect_true(all(startsWith(unlist(only), "glowspecific:")))
  expect_setequal(unname(c(without, unlist(only))), unname(all))
  expect_setequal(names(only), c("Disease-specific genes", "Class-specific genes"))
})

test_that("a user's own PanelApp download is used only by that user", {
  expect_identical(panelapp_dir("ana"), panels_dir)
  dir <- user_panelapp_dir("ana")
  dir.create(dir, recursive = TRUE)
  writeLines(c("panel_id\tpanel_name\tversion\tgene\tconfidence\tmode_of_inheritance",
               "1\tNew panel\t9.9\tGENEZ\t3\tBIALLELIC"), file.path(dir, "panelapp_genomics_england.tsv"))
  writeLines("2026-10-23", file.path(dir, "panelapp_date.txt"))
  expect_identical(panelapp_dir("ana"), dir)
  expect_identical(panelapp_dir("luis"), panels_dir)
  ana <- panel_catalog(panelapp_dir("ana"))
  expect_equal(ana$name[ana$id == "panelapp_genomics_england:1"], "New panel (v9.9)")
  expect_false("panelapp_genomics_england:1" %in% panel_catalog(panelapp_dir("luis"))$id &&
               identical(panel_catalog(panelapp_dir("luis"))$name[panel_catalog(panelapp_dir("luis"))$id == "panelapp_genomics_england:1"], "New panel (v9.9)"))
  expect_true(any(ana$source == "RetiGene"))   # RetiGene and the rest stay the lab ones
  expect_true(panelapp_info("ana")$own)
  expect_false(panelapp_info("luis")$own)
  unlink(file.path(shared_dir, "panelapp"), recursive = TRUE)
})

test_that("analysis rounds: saved per user and file, their variants can be hidden", {
  fresh_store()
  f <- tempfile(fileext = ".tsv"); writeLines("x", f)
  df <- make_variants(4)
  keys <- variant_key(df$CHROM, df$POS, df$REF, df$ALT)
  expect_length(list_rounds("ana", f), 0)
  add_round("ana", f, keys[1:2], "retina panel", "1 panel(s)")
  add_round("ana", f, keys[2:3])
  r <- list_rounds("ana", f)
  expect_equal(vapply(r, `[[`, "", "name"), c("retina panel", "Round 2"))
  expect_setequal(round_keys("ana", f), keys[1:3])
  expect_length(list_rounds("luis", f), 0)                         # per user
  expect_length(list_rounds("ana", tempfile()), 0)                 # per file
  expect_equal(filter_snv(df, snv_params(hidden_keys = round_keys("ana", f)))$POS, df$POS[4])
  delete_round("ana", f, 1)
  expect_equal(vapply(list_rounds("ana", f), `[[`, "", "name"), "Round 2")
  expect_error(add_round("ana", f, character(0)), "no variants")
})

test_that("default name of a round: sample and main filters", {
  expect_equal(round_name("P001", panels = c("Retinal disorders (v4.2)", "RP", "LCA"), max_af = 0.01, min_cadd = 20),
               "P001 · Retinal disorders (v4.2), RP +1 · AF ≤ 0.01 · CADD ≥ 20")
  expect_equal(round_name("P001"), "P001 · no filters")
  expect_equal(round_name("P001", hpo = 2, genes = c("A", "B")), "P001 \u00b7 2 HPO \u00b7 A, B")
})

test_that("PanelApp green, amber and red genes (confidence 3-4, 2, 0-1)", {
  dir <- tempfile(); dir.create(dir)
  writeLines(c("panel_id\tpanel_name\tversion\tgene\tconfidence\tmode_of_inheritance",
               "9\tDemo\t1.0\tG3\t3\tX", "9\tDemo\t1.0\tG4\t4\tX", "9\tDemo\t1.0\tA2\t2\tX", "9\tDemo\t1.0\tR1\t1\tX", "9\tDemo\t1.0\tR0\t0\tX"),
             file.path(dir, "panelapp_genomics_england.tsv"))
  cat <- panel_catalog(dir)
  id <- "panelapp_genomics_england:9"
  expect_setequal(panel_genes(id, catalog = cat), c("G3", "G4"))
  expect_setequal(panel_genes(id, amber = TRUE, catalog = cat), c("G3", "G4", "A2"))
  expect_setequal(panel_genes(id, red = TRUE, catalog = cat), c("G3", "G4", "A2", "R1", "R0"))
  expect_true(panelapp_has_red(cat))
})
