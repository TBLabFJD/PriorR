# Small helpers on variant tables.

# Variant key comparable between files: "chr" prefix optional, alleles in upper case
variant_ids <- function(chrom, pos, ref, alt) {
  paste(sub("^chr", "", as.character(chrom), ignore.case = TRUE), as.integer(pos),
        toupper(as.character(ref)), toupper(as.character(alt)), sep = "_")
}

# Key of each row as in the original VCF: from variant_id (chr_pos_ref_alt) when present,
# which keeps the VCF representation of indels, otherwise from CHROM/POS/REF/ALT
row_variant_ids <- function(DF) {
  if ("variant_id" %in% names(DF)) {
    parts <- strsplit(as.character(DF$variant_id), "_", fixed = TRUE)
    part <- function(i) vapply(parts, function(p) if (length(p) >= i) p[i] else NA_character_, "")
    variant_ids(part(1), part(2), part(3), part(4))
  } else {
    variant_ids(DF$CHROM, DF$POS, DF$REF, DF$ALT)
  }
}

# Picker filter where an empty value (NA) is kept when "NA" is selected
in_selection <- function(x, selected) {
  x %in% selected | (is.na(x) & any(is.na(selected) | selected == "NA"))
}

# Genes typed by the user: "BRCA1, BRCA2 TP53" -> c("BRCA1", "BRCA2", "TP53")
split_genes <- function(text) {
  genes <- trimws(unlist(strsplit(text %||% "", "[,;[:space:]]+")))
  unique(genes[genes != ""])
}

# For fields with several genes ("A,B;C"), does each row contain any of `genes`?
genes_overlap <- function(gene_field, genes) {
  parts <- strsplit(as.character(gene_field), "[,;|&]")
  vapply(parts, function(p) any(trimws(p) %in% genes), logical(1))
}

# Rows of DF inside any BED region (BED starts are 0-based: start < POS <= end)
in_bed <- function(DF, bed) {
  chrom <- sub("^chr", "", as.character(DF$CHROM), ignore.case = TRUE)
  bed_chrom <- sub("^chr", "", as.character(bed$chrom), ignore.case = TRUE)
  pos <- as.numeric(DF$POS)
  inside <- rep(FALSE, nrow(DF))
  for (ch in intersect(unique(chrom), unique(bed_chrom))) {
    rows <- which(chrom == ch)
    regions <- bed[bed_chrom == ch, ]
    regions <- regions[order(regions$start), ]
    # max end so far, to handle overlapping regions with findInterval
    max_end <- cummax(regions$end)
    i <- findInterval(pos[rows], regions$start, left.open = TRUE)
    inside[rows] <- i > 0 & pos[rows] <= max_end[pmax(i, 1)]
  }
  inside
}

`%||%` <- function(a, b) if (is.null(a)) b else a
