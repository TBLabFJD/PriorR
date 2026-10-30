# Cards of Variant information about the gene of the variant: the panels of the library that
# have it, and the other variants of the gene in this sample (compound heterozygotes).

# Panels of the library with the gene: source, name and level (green / amber / red for PanelApp,
# "in" for the other lists)
gene_panels <- function(gene, catalog = panel_catalog()) {
  empty <- data.frame(source = character(0), name = character(0), level = character(0))
  if (is.null(gene) || is.na(gene) || !nzchar(gene) || !nrow(catalog)) return(empty)
  has <- function(col) if (col %in% names(catalog)) vapply(catalog[[col]], function(g) gene %in% g, TRUE) else rep(FALSE, nrow(catalog))
  green <- has("genes"); amber <- has("amber"); red <- has("red")
  keep <- green | amber | red
  if (!any(keep)) return(empty)
  panelapp <- grepl("^PanelApp", catalog$source)
  level <- ifelse(red, "red", ifelse(amber, "amber", ifelse(panelapp, "green", "in")))
  out <- data.frame(source = catalog$source, name = catalog$name, level = level, stringsAsFactors = FALSE)[keep, , drop = FALSE]
  out[order(match(out$level, c("green", "in", "amber", "red")), out$source, out$name), , drop = FALSE]
}

# Genotype with an alternative allele (0/1, 1/1, 1|0, 1...)
has_alt <- function(gt) !is.na(gt) & grepl("[1-9]", gt)

# Parental origin of a variant of the proband of the family (PED / trio): "father", "mother",
# "both parents", "not in the parents" (de novo?) or NA (no family or genotypes)
variant_origin <- function(row, ped) {
  if (is.null(ped) || !nrow(ped)) return(NA_character_)
  p <- ped[ped$V6 == "2", , drop = FALSE]
  if (!nrow(p)) return(NA_character_)
  gt <- function(id) { col <- paste0(id, "_GT"); if (!is.na(id) && id != "0" && col %in% names(row)) as.character(row[[col]]) else NA }
  f <- gt(p$V3[1]); m <- gt(p$V4[1])
  if (is.na(f) && is.na(m)) return(NA_character_)
  fa <- has_alt(f); ma <- has_alt(m)
  if (fa && ma) "both parents" else if (fa) "father" else if (ma) "mother"
  else if (!is.na(f) && !is.na(m)) "not in the parents" else NA_character_
}

# The other variants of the gene of `v` in the table shown (one row per variant): position,
# change (HGVSp, else HGVSc, else consequence), genotypes of the samples and parental origin
gene_other_variants <- function(table, v, ped = NULL) {
  empty <- data.frame(position = character(0), change = character(0), genotypes = character(0), origin = character(0))
  if (is.null(table) || !nrow(table) || is.null(v$SYMBOL) || is.na(v$SYMBOL)) return(empty)
  same <- table[table$SYMBOL %in% v$SYMBOL, , drop = FALSE]
  key <- variant_key(same$CHROM, same$POS, same$REF, same$ALT)
  same <- same[!duplicated(key) & key != variant_key(v$CHROM, v$POS, v$REF, v$ALT), , drop = FALSE]
  if (!nrow(same)) return(empty)
  short <- function(x) sub("^[^:]*:", "", as.character(x))
  change <- if ("HGVSp" %in% names(same)) short(same$HGVSp) else rep(NA, nrow(same))
  if ("HGVSc" %in% names(same)) change[is.na(change) | !nzchar(change)] <- short(same$HGVSc)[is.na(change) | !nzchar(change)]
  if ("Consequence" %in% names(same)) change[is.na(change) | !nzchar(change)] <- gsub("_", " ", same$Consequence)[is.na(change) | !nzchar(change)]
  samples <- sample_names(names(same))
  genotypes <- if (length(samples)) apply(same[, paste0(samples, "_GT"), drop = FALSE], 1, function(g) paste(samples, g, collapse = " · ")) else ""
  data.frame(position = sprintf("%s:%s %s>%s", same$CHROM, key_text(same$POS), same$REF, same$ALT),
             change = change, genotypes = unname(genotypes),
             origin = vapply(seq_len(nrow(same)), function(i) variant_origin(same[i, , drop = FALSE], ped), ""),
             stringsAsFactors = FALSE)
}
