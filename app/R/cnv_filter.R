# CNV filtering, a pure function like filter_snv(). Each filter is applied only when its column
# is in the file (AnnotSV files differ between pipelines).

cnv_params <- function(...) {
  p <- list(samples = NULL, annotation_mode = NULL, type = NULL, region = NULL, min_acmg = NULL,
            min_programs = NULL, genes = character(0), panel_genes = NULL, excluded_genes = NULL,
            max_gain_af = NULL, max_loss_af = NULL, order_by_class = FALSE, columns = NULL,
            glow_sources = NULL, glow_max = NULL, glow_order = FALSE)   # GLOWgenes rankings (glowgenes.R)
  modifyList(p, list(...))
}

filter_cnv <- function(df, p = cnv_params()) {
  DF <- df
  has <- function(col) col %in% names(DF)
  num <- function(col) suppressWarnings(as.numeric(DF[[col]]))
  keep <- function(cond) DF[cond %in% TRUE, , drop = FALSE]

  if (!is.null(p$samples) && has("Samples_ID"))
    DF <- keep(vapply(strsplit(as.character(DF$Samples_ID), ",", fixed = TRUE), function(s) any(trimws(s) %in% p$samples), logical(1)))
  if (!is.null(p$annotation_mode) && has("Annotation_mode")) DF <- keep(DF$Annotation_mode %in% p$annotation_mode)
  if (!is.null(p$type) && has("SV_type")) DF <- keep(DF$SV_type %in% p$type)
  if (!is.null(p$region) && has("Location2")) DF <- keep(in_selection(DF$Location2, p$region))
  if (!is.null(p$min_acmg) && has("ACMG_class")) DF <- keep(num("ACMG_class") >= as.numeric(p$min_acmg) | is.na(DF$ACMG_class))
  if (!is.null(p$min_programs) && has("N_PROGRAMS")) DF <- keep(num("N_PROGRAMS") >= as.numeric(p$min_programs))
  # Gene_name can hold several genes: match whole gene names, not substrings
  if (has("Gene_name")) {
    if (length(p$genes)) DF <- keep(genes_overlap(DF$Gene_name, p$genes))
    if (!is.null(p$panel_genes)) DF <- keep(genes_overlap(DF$Gene_name, p$panel_genes))
    if (!is.null(p$excluded_genes)) DF <- keep(!genes_overlap(DF$Gene_name, p$excluded_genes))
  }
  if (!is.null(p$max_gain_af) && has("B_gain_AFmax")) DF <- keep(num("B_gain_AFmax") <= p$max_gain_af | is.na(DF$B_gain_AFmax))
  if (!is.null(p$max_loss_af) && has("B_loss_AFmax")) DF <- keep(num("B_loss_AFmax") <= p$max_loss_af | is.na(DF$B_loss_AFmax))

  # GLOWgenes: best rank of the genes of each CNV
  if (length(p$glow_sources) && has("Gene_name")) {
    DF$GLOW_rank <- cnv_glow_rank(DF$Gene_name, p$glow_sources)
    if (!is.null(p$glow_max)) DF <- keep(DF$GLOW_rank <= p$glow_max)
  }

  # Most pathogenic first: ACMG class, then AnnotSV ranking score (stable)
  if (isTRUE(p$order_by_class) && has("ACMG_class")) {
    score <- if (has("AnnotSV_ranking_score")) num("AnnotSV_ranking_score") else 0
    DF <- DF[order(-num("ACMG_class"), -score, na.last = TRUE), , drop = FALSE]
  }
  if (isTRUE(p$glow_order) && has("GLOW_rank")) DF <- DF[order(DF$GLOW_rank, na.last = TRUE), , drop = FALSE]

  out <- as.data.frame(DF)
  rownames(out) <- NULL
  columns <- if (is.null(p$columns)) names(out) else intersect(p$columns, names(out))
  attr(out, "display_columns") <- unique(c(intersect("GLOW_rank", names(out)), columns))
  out
}
