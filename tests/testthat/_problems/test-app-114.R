# Extracted from test-app.R:114

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
save_custom_panel("Mi panel", c("GENE2", "GENE3"))
me <- reactiveVal({ create_user("ana", "secret123"); authenticate("ana", "secret123") })
df <- make_variants(4, SYMBOL = c("GENE1", "GENE2", "GENE3", "GENE4"))
testServer(snv_server, args = list(user = me), {
    session$setInputs(consequence = NULL, region = NULL, type = NULL, canonical = FALSE, exclude_predicted = FALSE,
                      genes = "", max_gnomadg_af = 1, max_gnomadg_popmax = 1, max_gnomade_af = 1, max_fjd_af = 1,
                      min_cadd = "0", exclude_benign = FALSE, inheritance = "", glow = -1, orphanet = NULL,
                      panels = "custom:Mi panel", amber = FALSE, excluded_panels = NULL, phenotype = NULL, proband = "",
                      other_test_filter = FALSE, prioritize = FALSE, annotate = FALSE, columns = NULL)
    loader$loaded(list(df = df, profile = builtin_profiles()$None, restoring = TRUE))
    expect_equal(filtered()$SYMBOL, c("GENE2", "GENE3"))
  })
