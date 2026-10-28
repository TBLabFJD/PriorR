# Pathogenicity score filters: only the scores that the loaded file has (or that PriorR adds, like
# AlphaMissense and the gnomAD constraint) are shown. A missing score never removes a row when all
# the filters must pass; with "at least one", a row needs one of the chosen scores passing.

score_filter_defs <- data.frame(
  col = c("REVELScore", "AlphaMissense", "MutScore", "SpliceAI_DS_Max", "ada_score", "rf_score",
          "N_Pathogenic_pred", "Domino_Score", "phyloP30way_mammalian", "gnomAD_LOEUF"),
  label = c("Min REVEL", "Min AlphaMissense", "Min MutScore", "Min SpliceAI (max delta score)", "Min dbscSNV ADA",
            "Min dbscSNV RF", "Min number of pathogenic predictors", "Min DOMINO (dominant)", "Min phyloP (mammals)",
            "Max gene LOEUF (gnomAD)"),
  op = c(">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", "<="),
  min = c(0, 0, 0, 0, 0, 0, 0, 0, -10, 0), max = c(1, 1, 1, 1, 1, 1, 30, 1, 10, 2),
  step = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 1, 0.05, 0.5, 0.05),
  stringsAsFactors = FALSE)

# Filters for the columns of a table (AlphaMissense: the column of the file or the one added by PriorR)
score_filters_for <- function(columns) {
  defs <- score_filter_defs
  am <- alphamissense_column(columns)
  if (!is.null(am)) defs$col[defs$col == "AlphaMissense"] <- am
  defs[defs$col %in% columns, , drop = FALSE]
}

# Value of a slider that doesn't filter (min for ">=", max for "<=")
score_off_value <- function(def) if (def$op == ">=") def$min else def$max

# Rows passing the active score filters: scores = list of list(col, op, value)
score_keep <- function(DF, scores, mode = "all") {
  if (!length(scores)) return(rep(TRUE, nrow(DF)))
  tests <- lapply(scores, function(s) {
    x <- if (s$col %in% names(DF)) suppressWarnings(as.numeric(DF[[s$col]])) else rep(NA_real_, nrow(DF))
    pass <- if (s$op == ">=") x >= s$value else x <= s$value
    list(pass = pass %in% TRUE, missing = is.na(x))
  })
  if (identical(mode, "any")) Reduce(`|`, lapply(tests, `[[`, "pass"))
  else Reduce(`&`, lapply(tests, function(t) t$pass | t$missing))
}
