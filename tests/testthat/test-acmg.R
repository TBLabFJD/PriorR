test_that("ACMG suggestions from the annotation", {
  crit <- function(...) acmg_suggest(make_variants(1, ...)[1, ])$criterion
  # stop gained, absent from gnomAD -> PVS1 + PM2_Supporting = 9 points, likely pathogenic
  stop <- crit(Consequence = "stop_gained", gnomADg_AF_popmax = NA, gnomADe_AF_popmax = NA, gnomADg_AF = NA, gnomADe_AF = NA)
  expect_setequal(stop, c("PVS1", "PM2_Supporting"))
  expect_equal(acmg_class(stop), "Likely pathogenic")
  # missense with REVEL 0.95 -> PP3_Strong
  expect_true("PP3_Strong" %in% crit(REVELScore = 0.95))
  expect_true("PP3_Moderate" %in% crit(REVELScore = 0.8))
  expect_true("BP4" %in% crit(REVELScore = 0.25, SpliceAI_DS_Max = 0.01))
  # no REVEL: CADD
  expect_true("PP3" %in% crit(CADD_PHRED = 26))
  # common: BA1 is benign alone
  common <- crit(gnomADg_AF_popmax = 0.2)
  expect_true("BA1" %in% common)
  expect_equal(acmg_class(common), "Benign")
  # synonymous without splicing effect -> BP7
  expect_true("BP7" %in% crit(Consequence = "synonymous_variant", SpliceAI_DS_Max = 0.02))
  # intronic with SpliceAI 0.5 -> PP3
  expect_true("PP3" %in% crit(Consequence = "intron_variant", SpliceAI_DS_Max = 0.5))
  # de novo in the trio
  expect_true("PM6" %in% acmg_suggest(make_variants(1)[1, ], de_novo = TRUE)$criterion)
})

test_that("ACMG points and classes (Tavtigian 2020)", {
  expect_equal(acmg_total(c("PVS1", "PM2_Supporting")), 9)
  expect_equal(acmg_class(c("PVS1", "PS1")), "Pathogenic")
  expect_equal(acmg_class(c("PM2_Supporting", "PP3")), "Uncertain significance")
  expect_equal(acmg_class(c("BS1", "BP4")), "Likely benign")
  expect_equal(acmg_class(c("BS1", "BS2")), "Benign")
  expect_equal(acmg_class(character(0)), "Uncertain significance")
})

test_that("reanalysis: classified variants whose ClinVar, frequency or ACMG suggestion changed", {
  fresh_store()
  create_user("ana", "secret123")
  df <- make_variants(3, CLNSIG = c(NA, "Benign", NA), gnomADg_AF_popmax = c(0.00001, 0.2, 0.00001), REVELScore = c(0.1, 0.1, 0.1))
  key <- function(i) variant_key(df$CHROM[i], df$POS[i], df$REF[i], df$ALT[i])
  add_comment("ana", key(1), gene = "GENE1", classification = "Likely benign", snapshot = variant_snapshot(df[1, ]))
  add_comment("ana", key(2), gene = "GENE1", classification = "Benign", snapshot = variant_snapshot(df[2, ]))
  # variant 3 classified as VUS; the annotation now suggests more
  add_comment("ana", key(3), gene = "GENE1", classification = "Uncertain significance", snapshot = variant_snapshot(df[3, ]))
  expect_equal(nrow(reanalysis(df, visible_comments("ana"))), 0)   # nothing changed
  # later: ClinVar of variant 1 says pathogenic; variant 3 is now a stop gained
  df$CLNSIG[1] <- "Pathogenic"
  df$Consequence[3] <- "stop_gained"
  r <- reanalysis(df, visible_comments("ana"))
  expect_setequal(r$variant_key, c(key(1), key(3)))
  expect_match(r$changes[r$variant_key == key(1)], "ClinVar: not in ClinVar → Pathogenic")
  expect_match(r$changes[r$variant_key == key(3)], "ACMG suggestion: Likely benign → Likely pathogenic")
  expect_equal(unlist(all_comments()$snapshot[[1]]$CLNSIG), NULL)   # no ClinVar when classified
})

test_that("PS1 / PM5: pathogenic missense of ClinVar at the same residue (same transcript)", {
  f <- tempfile(fileext = ".tsv.gz"); con <- gzfile(f, "w")
  utils::write.table(data.frame(transcript = "NM_000350", aa_pos = c(1961, 1961), ref_aa = "Gly", alt_aa = c("Glu", "Arg"),
                                gene = "ABCA4", chrom = "chr1", pos = c(93992834, 93992835), ref = "C", alt = c("T", "G"),
                                clnsig = "Pathogenic"), con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  cv <- clinvar_missense(f)
  row <- function(...) make_variants(1, SYMBOL = "ABCA4", Feature = "NM_000350.3", CHROM = "chr1", POS = 93992900,
                                     REF = "A", ALT = "G", ...)[1, ]
  # same change p.Gly1961Glu by another nucleotide -> PS1
  crit <- acmg_suggest(row(HGVSp = "NP_000341.2:p.Gly1961Glu"), clinvar = cv)$criterion
  expect_true("PS1" %in% crit); expect_false("PM5" %in% crit)
  # another change at the residue -> PM5
  expect_true("PM5" %in% acmg_suggest(row(HGVSp = "NP_000341.2:p.Gly1961Val"), clinvar = cv)$criterion)
  # other transcript or other residue -> nothing
  other_tx <- make_variants(1, Feature = "NM_999999.1", HGVSp = "NP_x:p.Gly1961Glu")[1, ]
  expect_false(any(c("PS1", "PM5") %in% acmg_suggest(other_tx, clinvar = cv)$criterion))
  expect_false(any(c("PS1", "PM5") %in% acmg_suggest(row(HGVSp = "NP_000341.2:p.Gly1962Glu"), clinvar = cv)$criterion))
  # the very same variant of ClinVar is not PS1 (that is its own ClinVar classification)
  same <- make_variants(1, Feature = "NM_000350.3", CHROM = "chr1", POS = 93992834, REF = "C", ALT = "T", HGVSp = "NP_000341.2:p.Gly1961Glu")[1, ]
  crit <- acmg_suggest(same, clinvar = cv)$criterion
  expect_false("PS1" %in% crit); expect_true("PM5" %in% crit)
  expect_equal(acmg_total(c("PS1", "PM5_Supporting")), 5)
})
