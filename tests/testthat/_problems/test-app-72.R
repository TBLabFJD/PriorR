# Extracted from test-app.R:72

# prequel ----------------------------------------------------------------------
wait_loaded <- function(session, loaded) {
  for (i in 1:60) {
    if (!is.null(loaded())) return(loaded())
    Sys.sleep(0.5)
    session$elapse(700)
  }
  stop("not loaded")
}

# test -------------------------------------------------------------------------
fresh_store()
create_user("ana", "secret123")
me <- reactiveVal(authenticate("ana", "secret123"))
df <- make_variants(6, CADD_PHRED = c(30, 5, 20, NA, 16, 40), CLNSIG = c(NA, NA, "Benign", NA, "Pathogenic", NA),
                      Genomic_region = c("EXONIC", "EXONIC", "EXONIC", "INTRONIC", "SPLICING", NA))
testServer(snv_server, args = list(user = me), {
    session$setInputs(consequence = c("missense_variant", "NA"), region = c(genomic_regions, "NA"),
                      type = "SNV", canonical = FALSE, exclude_predicted = FALSE, genes = "",
                      max_gnomadg_af = 1, max_gnomadg_popmax = 1, max_gnomade_af = 1, max_fjd_af = 1,
                      min_cadd = "0", exclude_benign = FALSE, inheritance = "", glow = -1, orphanet = NULL,
                      panels = NULL, amber = FALSE, excluded_panels = NULL, phenotype = NULL, proband = "",
                      other_test_filter = FALSE, prioritize = FALSE, annotate = FALSE, columns = NULL)
    loader$loaded(list(df = df, profile = builtin_profiles()$None, restoring = FALSE))
    expect_equal(nrow(filtered()), 6)                 # default: no filters
    expect_equal(applied_profile(), "None")

    # the WES profile as the sidebar sets it (CADD >= 15 as a number, benign out, coding/splicing)
    session$setInputs(region = c("EXONIC", "SPLICING"), min_cadd = "15", exclude_benign = TRUE, prioritize = TRUE,
                      columns = c("SYMBOL", "CADD_PHRED"))
    out <- filtered()
    expect_equal(out$CADD_PHRED, c(16, 30))
    expect_equal(attr(out, "display_columns"), c("CHROM", "POS", "REF", "ALT", "SYMBOL", "CADD_PHRED"))

    # the current filters and columns are saved as a profile of the user
    session$setInputs(profile_name = "Mi retina", save_profile = 1)
    saved <- user_profiles("ana")[["Mi retina"]]
    expect_equal(unlist(saved$filters$region), c("EXONIC", "SPLICING"))
    expect_equal(saved$filters$min_cadd, 15)
    expect_equal(unlist(saved$columns), c("SYMBOL", "CADD_PHRED"))
    expect_null(saved$filters$type)                   # all types selected = no filter

    # Classify row 2 -> the comment is saved for that variant
    session$setInputs(classify = 2)
    info$open(list(key = variant_key(out$CHROM[2], out$POS[2], out$REF[2], out$ALT[2]), gene = "GENE1", label = "x"))
    session$setInputs(`info-comments-classification` = "Likely pathogenic", `info-comments-comment` = "check",
                      `info-comments-share` = "groups", `info-comments-save` = 1)
    session$setInputs(annotate = TRUE)
    out <- filtered()
    expect_equal(out$Classification, c(NA, "Likely pathogenic"))
    expect_equal(out$Comments, c(0L, 1L))
  })
