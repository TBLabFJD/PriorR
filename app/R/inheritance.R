# Inheritance patterns and ROH.

# Genotype tests (unphased or phased; haploid calls for chrX in males)
gt_het <- function(x) x %in% c("0/1", "0|1", "1|0")
gt_ref <- function(x) x %in% c("0/0", "0|0", "0")
gt_hom_alt <- function(x) x %in% c("1/1", "1|1", "1")
gt_carrier <- function(x) gt_het(x) | gt_hom_alt(x)

# TRUE for the rows where every column in `cols` passes `test` (TRUE when there are no columns)
all_gt <- function(DF, cols, test) {
  ok <- rep(TRUE, nrow(DF))
  for (col in cols) ok <- ok & test(as.character(DF[[col]]))
  ok
}

# Samples of a variant file: "P001" for the columns P001_GT, P001_DP... The columns of
# each variant caller (P001_DV_GT, P001_GK_GT...) are not samples but callers of a sample.
sample_names <- function(columns) {
  prefixes <- sub("_GT$", "", grep("_GT$", columns, value = TRUE))
  prefixes[!vapply(prefixes, function(p) any(startsWith(p, paste0(setdiff(prefixes, p), "_"))), logical(1))]
}

# Callers of a sample: c("DV", "DR", "GK") for the columns P001_DV_GT...
sample_callers <- function(columns, sample) {
  prefix <- paste0(sample, "_")
  gt <- grep("_GT$", columns, value = TRUE)
  callers <- gt[startsWith(gt, prefix) & gt != paste0(sample, "_GT")]
  sub("_GT$", "", substring(callers, nchar(prefix) + 1))
}

# Genotype column of the sample (the first sample)
sample_gt_column <- function(columns) {
  paste0(sample_names(columns)[1], "_GT")
}

# Genotype columns of each role in a PED file (family, id, father, mother, sex 1=male 2=female,
# phenotype 2=affected 1=unaffected). Only samples with a *_GT column in the data are used.
ped_samples <- function(ped, columns) {
  ped <- data.frame(lapply(ped, as.character), stringsAsFactors = FALSE)
  gt_col <- function(ids) {
    cols <- paste0(ids, "_GT")
    cols[!is.na(ids) & ids != "0" & cols %in% columns]
  }
  affected <- ped$V6 == "2"
  unaffected <- ped$V6 == "1"
  male <- ped$V5 == "1"
  first_affected <- ped[affected, ][1, ]
  father <- gt_col(first_affected$V3)
  mother <- gt_col(first_affected$V4)
  list(
    affected = gt_col(ped$V2[affected]),
    unaffected = gt_col(ped$V2[unaffected]),
    father = if (length(father)) father else NULL,
    mother = if (length(mother)) mother else NULL,
    relatives = setdiff(gt_col(ped$V2[unaffected]), c(father, mother)),
    affected_males = gt_col(ped$V2[affected & male]),
    affected_females = gt_col(ped$V2[affected & !male]),
    unaffected_males = gt_col(ped$V2[unaffected & male])
  )
}

# Compound heterozygous candidates: genes with at least two different heterozygous
# variants in the affected samples (recessive OMIM genes). With both parents' genotypes,
# at least one variant has to come from the father and another one from the mother.
compound_het <- function(DF, gt_affected, gt_father = NULL, gt_mother = NULL, gt_unaffected = character(0)) {
  gt <- function(col) as.character(DF[[col]])

  keep <- all_gt(DF, gt_affected, gt_het) &
    all_gt(DF, gt_unaffected, function(x) !gt_hom_alt(x)) &
    grepl("recessive", DF$OMIM_phenotype, ignore.case = TRUE) &
    !is.na(DF$SYMBOL)

  variant <- row_variant_ids(DF)
  gene <- as.character(DF$SYMBOL)

  if (!is.null(gt_father) && !is.null(gt_mother)) {
    paternal <- keep & gt_het(gt(gt_father)) & gt_ref(gt(gt_mother))
    maternal <- keep & gt_het(gt(gt_mother)) & gt_ref(gt(gt_father))
    genes <- intersect(gene[paternal], gene[maternal])
    rows <- (paternal | maternal) & gene %in% genes
  } else {
    n_variants <- tapply(variant[keep], gene[keep], function(v) length(unique(v)))
    genes <- names(n_variants)[n_variants >= 2]
    rows <- keep & gene %in% genes
  }
  DF[rows, , drop = FALSE]
}

# Read depth and alternative reads of a sample (columns <id>_DP and <id>_AD "ref,alt" or "ref_alt");
# NA when the file doesn't have them
sample_dp <- function(DF, gt_col) {
  col <- sub("_GT$", "_DP", gt_col)
  if (col %in% names(DF)) suppressWarnings(as.numeric(DF[[col]])) else rep(NA_real_, nrow(DF))
}
sample_alt_reads <- function(DF, gt_col) {
  col <- sub("_GT$", "_AD", gt_col)
  if (!col %in% names(DF)) return(rep(NA_real_, nrow(DF)))
  ad <- strsplit(as.character(DF[[col]]), "[,_]")
  vapply(ad, function(x) if (length(x) >= 2) sum(suppressWarnings(as.numeric(x[-1]))) else NA_real_, 0)
}
sample_vaf <- function(DF, gt_col) {
  col <- sub("_GT$", "_VAF", gt_col)
  if (col %in% names(DF)) return(suppressWarnings(as.numeric(DF[[col]])))
  alt <- sample_alt_reads(DF, gt_col)
  alt / sample_dp(DF, gt_col)
}

# De novo quality thresholds (only applied when the file has DP / AD / VAF)
denovo_settings <- list(min_dp = 10, max_parent_alt_reads = 1, proband_vaf = c(0.2, 0.8))

# De novo candidates of a trio: the proband carries the variant and both parents are homozygous
# reference with enough depth and (almost) no alternative reads. chrX in a male proband: only the
# mother is checked; chrY: only the father. Parents without a call ("./.") never give a de novo.
de_novo <- function(DF, proband, father, mother, proband_male = FALSE, s = denovo_settings) {
  gt <- function(col) as.character(DF[[col]])
  chrom <- sub("^chr", "", DF$CHROM)
  on_x <- chrom == "X"
  on_y <- chrom == "Y"
  ok_dp <- function(col) { dp <- sample_dp(DF, col); is.na(dp) | dp >= s$min_dp }
  clean_parent <- function(col) gt_ref(gt(col)) & ok_dp(col) &
    (is.na(sample_alt_reads(DF, col)) | sample_alt_reads(DF, col) <= s$max_parent_alt_reads)
  vaf <- sample_vaf(DF, proband)
  het_ok <- gt_het(gt(proband)) & (is.na(vaf) | (vaf >= s$proband_vaf[1] & vaf <= s$proband_vaf[2]))
  carrier <- ok_dp(proband) & (het_ok | (proband_male & (on_x | on_y) & gt_hom_alt(gt(proband))))
  parents <- ifelse(proband_male & on_x, clean_parent(mother),
             ifelse(proband_male & on_y, clean_parent(father), clean_parent(father) & clean_parent(mother)))
  carrier & parents
}

# A trio as a PED table (V1 family ... V6 phenotype) from the samples chosen in the sidebar
trio_ped <- function(proband, father, mother, proband_sex = "1") {
  ids <- c(proband, father, mother)
  ids[is.na(ids) | !nzchar(ids)] <- "0"
  rows <- data.frame(V1 = "trio", V2 = ids, V3 = c(ids[2], "0", "0"), V4 = c(ids[3], "0", "0"),
                     V5 = c(proband_sex, "1", "2"), V6 = c("2", "1", "1"), stringsAsFactors = FALSE)
  rows[rows$V2 != "0", , drop = FALSE]
}

# Rows compatible with an inheritance mode. Without a PED file the sample genotype column is used.
inheritance_filter <- function(DF, mode, ped = NULL) {
  if (is.null(mode) || mode == "") return(DF)
  dominant <- grepl("dominant", DF$OMIM_phenotype, ignore.case = TRUE)
  recessive <- grepl("recessive", DF$OMIM_phenotype, ignore.case = TRUE)
  on_x <- DF$CHROM %in% c("chrX", "X")

  if (is.null(ped)) {
    gt <- sample_gt_column(names(DF))
    g <- as.character(DF[[gt]])
    keep <- switch(mode,
      "De novo" = rep(FALSE, nrow(DF)),   # needs the parents (PED or trio)
      "Dominant" = gt_het(g) & dominant,
      "Recessive-Homozygous" = (gt_hom_alt(g) | g == "1/2") & recessive,
      "Recessive-Compound heterozygous" = return(compound_het(DF, gt)),
      "X-linked" = on_x)
    return(DF[keep %in% TRUE, , drop = FALSE])
  }

  fam <- ped_samples(ped, names(DF))
  keep <- switch(mode,
    "De novo" = if (length(fam$affected) && !is.null(fam$father) && !is.null(fam$mother)) {
      proband <- fam$affected[1]
      male <- ped$V5[match(sub("_GT$", "", proband), ped$V2)] == "1"
      de_novo(DF, proband, fam$father, fam$mother, proband_male = isTRUE(male))
    } else rep(FALSE, nrow(DF)),
    # every affected heterozygous, every unaffected homozygous reference
    "Dominant" = all_gt(DF, fam$affected, gt_het) & all_gt(DF, fam$unaffected, gt_ref) & dominant,
    # every affected homozygous alternative, parents (when available) carriers,
    # other unaffected relatives not homozygous alternative
    "Recessive-Homozygous" = all_gt(DF, fam$affected, gt_hom_alt) &
      all_gt(DF, c(fam$father, fam$mother), gt_het) &
      all_gt(DF, fam$relatives, function(x) !gt_hom_alt(x)) & recessive,
    "Recessive-Compound heterozygous" =
      return(compound_het(DF, fam$affected, fam$father, fam$mother, fam$relatives)),
    # affected males hemizygous, affected females carriers, unaffected males without the variant
    "X-linked" = on_x & all_gt(DF, fam$affected_males, gt_hom_alt) &
      all_gt(DF, fam$affected_females, gt_carrier) &
      all_gt(DF, fam$unaffected_males, function(x) !gt_carrier(x)))
  DF[keep %in% TRUE, , drop = FALSE]
}

# Rows inside a ROH of the given samples (all of them); with no samples, the only _ROH column
roh_rows <- function(DF, samples = NULL) {
  if (length(samples) == 0) {
    samples <- sub("_ROH$", "", names(DF)[endsWith(names(DF), "_ROH")])
    if (length(samples) != 1) return(rep(FALSE, nrow(DF)))
  }
  cols <- intersect(paste0(samples, "_ROH"), names(DF))
  if (length(cols) == 0) return(rep(FALSE, nrow(DF)))
  rowSums(as.matrix(DF[, cols, drop = FALSE]) == "True", na.rm = TRUE) == length(cols)
}
