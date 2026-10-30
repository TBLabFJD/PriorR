# SNV filtering: a pure function (data + parameters -> table), so the same parameters
# always give the same table and it can be tested without Shiny.

# Default parameters (nothing filtered)
snv_params <- function(...) {
  p <- list(
    consequence = NULL, region = NULL, type = NULL,   # NULL = no filter
    canonical = FALSE, exclude_predicted = FALSE, genes = character(0),
    max_af = list(), min_cadd = NULL,
    scores = NULL, scores_mode = "all",   # score filters (score_filters.R)
    inheritance = "", ped = NULL, glow_max = NULL, glow_sources = NULL,   # max GLOWgenes rank (column or rankings)
    panel_genes = NULL, excluded_genes = NULL, orphanet = NULL, bed = NULL,
    phenotype_genes = NULL, other_test_ids = NULL, hide_other_test = FALSE,
    exclude_clinvar_benign = FALSE, prioritize = FALSE,
    annotations = NULL, columns = NULL, roh_samples = NULL,
    only_keys = NULL,   # variant keys to keep (reanalysis)
    hidden_keys = NULL,   # variant keys to hide (previous analysis rounds, rounds.R)
    glow = NULL)   # list(sources): the GLOWgenes ranking view instead of the order above
  modifyList(p, list(...))
}

# Key of a variant as written in the table (used for comments and classifications)
variant_key <- function(chrom, pos, ref, alt) {
  na_empty <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  paste(na_empty(chrom), na_empty(pos), na_empty(ref), na_empty(alt), sep = ":")
}

filter_snv <- function(df, p = snv_params()) {
  DF <- df
  num <- function(col) suppressWarnings(as.numeric(DF[[col]]))
  keep <- function(cond) DF[cond %in% TRUE, , drop = FALSE]

  if (!is.null(p$consequence)) DF <- keep(in_selection(DF$Consequence, p$consequence))
  if (!is.null(p$region)) DF <- keep(in_selection(DF$Genomic_region, p$region))
  if (!is.null(p$type)) DF <- keep(DF$VARIANT_CLASS %in% p$type)
  if (isTRUE(p$canonical)) DF <- keep(DF$CANONICAL == "YES")
  # RefSeq predicted transcripts (XM_, XR_)
  if (isTRUE(p$exclude_predicted)) DF <- keep(!grepl("^X[MR]_", DF$Feature))
  if (length(p$genes)) DF <- keep(DF$SYMBOL %in% p$genes)
  for (col in names(p$max_af)) {
    # not a number ("", ".") = no frequency, kept like an empty value
    if (!is.null(p$max_af[[col]]) && col %in% names(DF)) DF <- keep(num(col) <= p$max_af[[col]] | is.na(num(col)))
  }
  if (!is.null(p$min_cadd)) DF <- keep(num("CADD_PHRED") >= as.numeric(p$min_cadd) | is.na(num("CADD_PHRED")))
  if (length(p$scores)) DF <- keep(score_keep(DF, p$scores, p$scores_mode))

  DF <- inheritance_filter(DF, p$inheritance, p$ped)

  if (!is.null(p$glow_max) && is.null(p$glow)) DF <- keep(glow_rank(DF, p$glow_sources) <= p$glow_max)
  if (!is.null(p$panel_genes)) DF <- keep(DF$SYMBOL %in% p$panel_genes)
  if (!is.null(p$excluded_genes)) DF <- keep(!DF$SYMBOL %in% p$excluded_genes)

  # Orphanet panel: order by proximity to the panel (GLOWgenes rank)
  if (length(p$orphanet)) {
    DF <- dplyr::left_join(DF, orphanet_rank(p$orphanet), by = "SYMBOL")
    DF <- DF[order(DF$Orphanet_GLOW_rank, na.last = TRUE), , drop = FALSE]
  }
  if (!is.null(p$bed)) DF <- keep(in_bed(DF, p$bed))

  # Classifications and comments visible to the user
  if (!is.null(p$annotations)) {
    key <- variant_key(DF$CHROM, DF$POS, DF$REF, DF$ALT)
    i <- match(key, p$annotations$variant_key)
    DF$Classification <- p$annotations$classification[i]
    DF$Comments <- ifelse(is.na(i), 0L, p$annotations$n_comments[i])
  }
  if (!is.null(p$phenotype_genes)) DF <- keep(DF$SYMBOL %in% p$phenotype_genes)
  if (!is.null(p$only_keys)) DF <- keep(variant_key(DF$CHROM, DF$POS, DF$REF, DF$ALT) %in% p$only_keys)
  if (length(p$hidden_keys)) DF <- keep(!variant_key(DF$CHROM, DF$POS, DF$REF, DF$ALT) %in% p$hidden_keys)

  in_other <- if (is.null(p$other_test_ids)) rep(FALSE, nrow(DF)) else row_variant_ids(DF) %in% p$other_test_ids
  if (isTRUE(p$hide_other_test)) {
    DF <- DF[!in_other, , drop = FALSE]
    in_other <- in_other[!in_other]
  }
  if (isTRUE(p$exclude_clinvar_benign)) {
    benign <- grepl("benign", DF$CLNSIG, ignore.case = TRUE)
    DF <- DF[!benign, , drop = FALSE]
    in_other <- in_other[!benign]
  }

  if (isTRUE(p$prioritize)) {
    o <- order(match(priority_class(DF), c("red", "orange"), nomatch = 3L))
    DF <- DF[o, , drop = FALSE]
    in_other <- in_other[o]
  }

  # GLOWgenes ranking view (glowgenes.R): its own pathogenicity filters and order
  if (!is.null(p$glow)) {
    keep_key <- row_variant_ids(DF)
    DF <- glow_view(DF, p$glow$sources)
    if (!is.null(p$glow_max)) DF <- DF[(DF$GLOW_rank <= p$glow_max) %in% TRUE, , drop = FALSE]
    # only new candidates: without the genes of the disease panel (rank 0)
    if (isTRUE(p$glow$hide_rank0)) DF <- DF[!(DF$GLOW_rank %in% 0), , drop = FALSE]
    in_other <- in_other[match(row_variant_ids(DF), keep_key)]
  }

  out <- as.data.frame(DF)
  rownames(out) <- NULL

  # Displayed columns (all the columns are kept for the variant details and the colors):
  # hiding a column never changes the filters. Any column can be hidden (the variant of a row
  # is taken from the whole table); the columns added above are always shown, and the variant
  # key when nothing else would be
  added <- setdiff(names(out), names(df))
  columns <- if (is.null(p$columns)) names(df) else intersect(p$columns, names(out))
  shown <- unique(c(intersect("GLOW_rank", names(out)), columns, added))
  if (!length(shown)) shown <- intersect(c("CHROM", "POS", "REF", "ALT"), names(out))
  attr(out, "display_columns") <- shown

  # Row colors
  attr(out, "priority") <- priority_class(out)
  attr(out, "roh") <- roh_rows(out, p$roh_samples)
  attr(out, "not_in_other_test") <- if (is.null(p$other_test_ids)) rep(FALSE, nrow(out)) else !in_other
  out
}

# The displayed part of a filtered table
display_table <- function(DF) {
  DF[, attr(DF, "display_columns"), drop = FALSE]
}
