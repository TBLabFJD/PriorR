# Variant prioritization.

# Prioritization class of each row: "red", "orange" or NA.
# Single source of truth for both the PRIORITIZATION ordering and the COLORING,
# so they always agree. A row that matches red and orange is red (never duplicated).
priority_class <- function(DF) {
  num <- function(col) suppressWarnings(as.numeric(as.character(DF[[col]])))
  is_true <- function(x) x %in% TRUE
  popmax <- num("gnomADg_AF_popmax")
  patho_ratio <- num("N_Pathogenic_pred") / num("N_predictions")
  n_pred <- num("N_predictions")
  cadd <- num("CADD_PHRED")
  mutscore <- num("MutScore")
  ada <- num("ada_score")
  rf <- num("rf_score")
  clnsig <- as.character(DF[["CLNSIG"]])
  consequence <- as.character(DF[["Consequence"]])
  splice_or_intron <- grepl("splice", consequence) | grepl("intron", consequence)
  patho_conseq <- c("start_lost", "stop_gained", "frameshift_variant", "stop_lost", "inframe_deletion", "inframe_insertion")
  clinvar_plp <- c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")
  # ClinVar renamed "interpretations" to "classifications" in 2023: treat both the same
  clinvar_conflicting <- c("Conflicting_interpretations_of_pathogenicity", "Conflicting_classifications_of_pathogenicity")
  rare <- popmax < 0.001 | is.na(popmax)
  low_freq <- popmax < 0.01 & popmax > 0.001

  red <- is_true(
    (
      clnsig %in% clinvar_plp |
        (rare & patho_ratio > 0.60 & cadd > 25 & (mutscore > 0.80 | is.na(mutscore))) |
        (rare & patho_ratio > 0.60 & (cadd > 25 | is.na(cadd)) & mutscore > 0.80) |
        (consequence %in% patho_conseq & rare & patho_ratio > 0.75 & is.na(cadd) & is.na(mutscore)) |
        (splice_or_intron & rare & (ada > 0.8 | rf > 0.8))
    ) &
      !grepl("benign", clnsig, ignore.case = TRUE) &
      !grepl("Uncertain_significance", clnsig, ignore.case = TRUE)
  )

  orange <- is_true(
    (low_freq & (patho_ratio > 0.60 | n_pred == 0) & cadd > 25 & (mutscore > 0.80 | is.na(mutscore))) |
      (low_freq & patho_ratio > 0.60 & (cadd > 25 | mutscore > 0.80)) |
      (splice_or_intron & low_freq & patho_ratio > 0.60 & (ada > 0.8 | rf > 0.8)) |
      (grepl("pathogenic", clnsig, ignore.case = TRUE) & !(clnsig %in% c(clinvar_conflicting, clinvar_plp))) |
      (rare & patho_ratio > 0.60 & (cadd > 25 | mutscore > 0.80) & clnsig %in% "Uncertain_significance")
  )

  ifelse(red, "red", ifelse(orange, "orange", NA_character_))
}

# Stable ordering: red rows, then orange, then the rest (each group keeps its order)
prioritize <- function(DF) {
  DF[order(match(priority_class(DF), c("red", "orange"), nomatch = 3L)), , drop = FALSE]
}

# Orphanet panels available ("retinal_dystrophy_panel_glow", ...)
orphanet_panels <- function(pattern = "_panel_glow.txt$") {
  tools::file_path_sans_ext(list.files(orphanet_dir, pattern = pattern))
}

# Gene ranking of the selected Orphanet panels (panel genes first, then the rest
# of the genome by GLOWgenes proximity). With several panels, the best rank wins.
orphanet_rank <- function(panels) {
  ranks <- lapply(panels, function(p) {
    f <- file.path(orphanet_dir, paste0(p, "_new.txt"))
    if (!file.exists(f)) f <- file.path(orphanet_dir, paste0(p, ".txt"))
    r <- data.table::fread(f, header = FALSE, sep = "\t", data.table = FALSE)
    data.frame(SYMBOL = as.character(r$V1), Orphanet_GLOW_rank = as.numeric(r$V3), stringsAsFactors = FALSE)
  })
  ranks <- do.call(rbind, ranks)
  stats::aggregate(Orphanet_GLOW_rank ~ SYMBOL, data = ranks, FUN = min)
}

# Genes of Orphanet panels (plain "_panel.txt" lists, used by CNVs)
orphanet_genes <- function(panels) {
  unique(unlist(lapply(panels, function(p) read_gene_list(file.path(orphanet_dir, paste0(p, ".txt"))))))
}

# What red and orange mean (shown above the table in the PriorR priority view and in the Help)
priority_rules <- list(
  red = c("ClinVar pathogenic or likely pathogenic",
          "rare (gnomAD genomes popmax < 0.001 or absent), > 60% of the predictors pathogenic and CADD > 25 and / or MutScore > 0.8",
          "rare start / stop lost, stop gained, frameshift or in-frame indel with > 75% pathogenic predictors (when there is no CADD or MutScore)",
          "rare splicing or intronic variant with dbscSNV ADA or RF > 0.8"),
  red_never = "never when ClinVar says benign, likely benign or uncertain significance",
  orange = c("low frequency (popmax 0.001-0.01), > 60% of the predictors pathogenic and CADD > 25 or MutScore > 0.8",
             "low frequency splicing or intronic variant with > 60% pathogenic predictors and ADA or RF > 0.8",
             "other pathogenic ClinVar annotations (e.g. low penetrance, risk factor), not conflicting",
             "rare ClinVar uncertain significance with > 60% pathogenic predictors and CADD > 25 or MutScore > 0.8"))
