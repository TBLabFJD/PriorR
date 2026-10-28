# Suggested ACMG/AMP criteria of a variant (Richards et al. 2015 with the ClinGen SVI updates), from
# the annotation of the file. They are only suggestions: the analyst confirms or removes them.
# Points (Tavtigian et al. 2020): very strong 8, strong 4, moderate 2, supporting 1 (benign negative).
#   >= 10 pathogenic, 6-9 likely pathogenic, 0-5 uncertain, -1 to -6 likely benign, <= -7 benign
# - PVS1: null variant (stop, frameshift, canonical splice site, start lost); check that loss of
#   function is the disease mechanism, the last exon / NMD and alternative transcripts
# - PM2_Supporting: absent or very rare in gnomAD (popmax <= 0.0001)
# - PM4: in-frame indel or stop lost
# - PP3 / BP4: REVEL (Pejaver et al. 2022) for missense, or CADD when there is no REVEL;
#   SpliceAI >= 0.2 (PP3) / <= 0.1 (BP4) for splicing (Walker et al. 2023)
# - BP7: synonymous with SpliceAI <= 0.1
# - BA1 / BS1: popmax > 0.05 / > 0.01 (BS1 depends on the disease)
# - PS1 / PM5: same amino acid change / other missense at the same residue, pathogenic in ClinVar
#   (same transcript; Program_files/ClinVar)
# - PM6: de novo in the trio (PS2 when the parentage is confirmed)

acmg_points <- c(VeryStrong = 8, Strong = 4, Moderate = 2, Supporting = 1)

acmg_strength_name <- function(code, strength) {
  default <- c(PVS1 = "VeryStrong", PS1 = "Strong", PM2 = "Moderate", PM4 = "Moderate", PM5 = "Moderate", PM6 = "Moderate", PP3 = "Supporting",
               BP4 = "Supporting", BP7 = "Supporting", BA1 = "StandAlone", BS1 = "Strong")
  if (identical(default[[code]], strength)) code else paste(code, strength, sep = "_")
}

# Pathogenic / likely pathogenic missense variants of ClinVar by transcript and residue
# (Program_files/ClinVar, scripts/update_clinvar_missense.R), for PS1 and PM5
.clinvar_cache <- new.env()
clinvar_missense <- function(path = file.path(data_dir, "ClinVar", "clinvar_missense_plp.tsv.gz")) {
  if (is.null(.clinvar_cache[[path]])) {
    t <- if (file.exists(path)) utils::read.delim(gzfile(path), quote = "", stringsAsFactors = FALSE) else data.frame()
    .clinvar_cache[[path]] <- list(table = t, index = if (nrow(t)) split(seq_len(nrow(t)), paste(t$transcript, t$aa_pos)) else list())
  }
  .clinvar_cache[[path]]
}

# PS1 / PM5 of a missense variant: same amino acid change (by another nucleotide change) or another
# change at the same residue, pathogenic in ClinVar, in the same transcript (NM_ without version)
clinvar_same_residue <- function(v, clinvar = clinvar_missense()) {
  hgvsp <- as.character(v$HGVSp %||% NA); feature <- as.character(v$Feature %||% NA)
  if (is.na(hgvsp) || is.na(feature) || !nrow(clinvar$table)) return(NULL)
  m <- regmatches(hgvsp, regexec("p\\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})$", hgvsp))[[1]]
  if (length(m) != 4) return(NULL)
  rows <- clinvar$index[[paste(sub("\\..*$", "", feature), m[3])]]
  if (!length(rows)) return(NULL)
  hits <- clinvar$table[rows, ]
  same_nt <- hits$chrom == sub("^(chr)?", "chr", as.character(v$CHROM)) & hits$pos == as.numeric(v$POS) &
    hits$ref == as.character(v$REF) & hits$alt == as.character(v$ALT)
  hits <- hits[!same_nt, ]                       # the variant itself is CLNSIG, not PS1
  label <- function(h) paste(sprintf("p.%s%d%s (%s)", h$ref_aa, h$aa_pos, h$alt_aa, h$clnsig), collapse = "; ")
  same_aa <- hits[hits$alt_aa == m[4], ]
  other_aa <- hits[hits$alt_aa != m[4], ]
  list(ps1 = if (nrow(same_aa)) label(same_aa[1, ]), pm5 = if (nrow(other_aa)) label(utils::head(unique(other_aa[, c("ref_aa", "aa_pos", "alt_aa", "clnsig")]), 3)))
}

# Criteria of one variant (a row of the table): data.frame(code, criterion, strength, points, reason)
acmg_suggest <- function(v, de_novo = FALSE, clinvar = clinvar_missense()) {
  num <- function(col) if (col %in% names(v)) suppressWarnings(as.numeric(v[[col]])) else NA_real_
  txt <- function(col) if (col %in% names(v)) as.character(v[[col]]) else NA_character_
  out <- list()
  add <- function(code, strength, reason) {
    benign <- startsWith(code, "B")
    pts <- if (strength == "StandAlone") -8 else acmg_points[[strength]] * if (benign) -1 else 1
    out[[length(out) + 1]] <<- data.frame(code = code, criterion = acmg_strength_name(code, strength), strength = strength,
                                          points = pts, reason = reason, stringsAsFactors = FALSE)
  }
  consequence <- tolower(txt("Consequence") %||% "")
  consequence[is.na(consequence)] <- ""
  popmax <- suppressWarnings(max(num("gnomADg_AF_popmax"), num("gnomADe_AF_popmax"), num("gnomADg_AF"), num("gnomADe_AF"), na.rm = TRUE))
  if (!is.finite(popmax)) popmax <- NA
  revel <- num("REVELScore"); cadd <- num("CADD_PHRED"); spliceai <- num("SpliceAI_DS_Max")
  missense <- grepl("missense", consequence)
  null_variant <- grepl("stop_gained|frameshift|splice_acceptor_variant|splice_donor_variant|start_lost|transcript_ablation", consequence)

  # Population
  if (!is.na(popmax) && popmax > 0.05) add("BA1", "StandAlone", sprintf("gnomAD popmax %.3g > 0.05", popmax))
  else if (!is.na(popmax) && popmax > 0.01) add("BS1", "Strong", sprintf("gnomAD popmax %.3g > 0.01 (check the disease prevalence)", popmax))
  else if (is.na(popmax) || popmax <= 1e-4) add("PM2", "Supporting", if (is.na(popmax)) "absent from gnomAD" else sprintf("gnomAD popmax %.3g <= 0.0001", popmax))

  # Variant type
  if (null_variant) {
    loeuf <- num("gnomAD_LOEUF")
    add("PVS1", "VeryStrong", paste0("null variant (", sub(",.*", "", consequence), ")",
                                     if (!is.na(loeuf)) sprintf("; gene LOEUF %.2f", loeuf),
                                     "; check LoF mechanism, last exon/NMD and transcripts"))
  }
  if (grepl("inframe_insertion|inframe_deletion|stop_lost", consequence) && !null_variant)
    add("PM4", "Moderate", "in-frame indel or stop lost (check it is not a repeat region)")

  # Computational evidence (only one of PP3 / BP4)
  if (missense && !is.na(revel)) {
    if (revel >= 0.932) add("PP3", "Strong", sprintf("REVEL %.3f >= 0.932", revel))
    else if (revel >= 0.773) add("PP3", "Moderate", sprintf("REVEL %.3f >= 0.773", revel))
    else if (revel >= 0.644) add("PP3", "Supporting", sprintf("REVEL %.3f >= 0.644", revel))
    else if (revel <= 0.052 && (is.na(spliceai) || spliceai <= 0.1)) add("BP4", "Strong", sprintf("REVEL %.3f <= 0.052", revel))
    else if (revel <= 0.183 && (is.na(spliceai) || spliceai <= 0.1)) add("BP4", "Moderate", sprintf("REVEL %.3f <= 0.183", revel))
    else if (revel <= 0.290 && (is.na(spliceai) || spliceai <= 0.1)) add("BP4", "Supporting", sprintf("REVEL %.3f <= 0.290", revel))
  } else if (missense && !is.na(cadd)) {
    if (cadd >= 28.1) add("PP3", "Moderate", sprintf("CADD %.1f >= 28.1 (no REVEL)", cadd))
    else if (cadd >= 25.3) add("PP3", "Supporting", sprintf("CADD %.1f >= 25.3 (no REVEL)", cadd))
    else if (cadd <= 17.3 && (is.na(spliceai) || spliceai <= 0.1)) add("BP4", "Moderate", sprintf("CADD %.1f <= 17.3 (no REVEL)", cadd))
    else if (cadd <= 22.7 && (is.na(spliceai) || spliceai <= 0.1)) add("BP4", "Supporting", sprintf("CADD %.1f <= 22.7 (no REVEL)", cadd))
  }
  has_pp3 <- any(vapply(out, function(x) x$code == "PP3", logical(1)))
  has_bp4 <- any(vapply(out, function(x) x$code == "BP4", logical(1)))
  if (!null_variant && !is.na(spliceai)) {
    if (spliceai >= 0.2 && !has_pp3) add("PP3", "Supporting", sprintf("SpliceAI %.2f >= 0.2", spliceai))
    else if (spliceai <= 0.1 && !missense && !has_bp4) {
      if (grepl("synonymous", consequence)) add("BP7", "Supporting", sprintf("synonymous, SpliceAI %.2f <= 0.1", spliceai))
      else if (grepl("intron|utr|upstream|downstream", consequence)) add("BP4", "Supporting", sprintf("SpliceAI %.2f <= 0.1", spliceai))
    }
  }

  # Same residue pathogenic in ClinVar (PS1 and PM5 are not used together)
  if (missense) {
    cv <- clinvar_same_residue(v, clinvar)
    if (!is.null(cv$ps1)) add("PS1", "Strong", paste("same amino acid change pathogenic in ClinVar:", cv$ps1, "(check splicing)"))
    else if (!is.null(cv$pm5)) add("PM5", "Moderate", paste("other pathogenic missense at this residue in ClinVar:", cv$pm5))
  }

  if (isTRUE(de_novo)) add("PM6", "Moderate", "de novo in the trio (PS2 if the parentage is confirmed)")

  if (!length(out)) return(data.frame(code = character(0), criterion = character(0), strength = character(0),
                                      points = numeric(0), reason = character(0)))
  do.call(rbind, out)
}

# Points of chosen criteria ("PM2_Supporting", "PP3", ...); BA1 alone gives benign
acmg_total <- function(criteria) {
  if (!length(criteria)) return(0)
  if (any(startsWith(criteria, "BA1"))) return(-8)
  sum(vapply(criteria, function(cr) {
    parts <- strsplit(cr, "_")[[1]]
    strength <- if (length(parts) > 1) parts[2] else switch(substr(parts[1], 1, 2), PV = "VeryStrong", PS = "Strong",
                                                             PM = "Moderate", BS = "Strong", "Supporting")
    acmg_points[[strength]] * if (startsWith(parts[1], "B")) -1 else 1
  }, 0))
}

acmg_class <- function(criteria) {
  if (any(startsWith(criteria, "BA1"))) return("Benign")
  total <- acmg_total(criteria)
  if (total >= 10) "Pathogenic" else if (total >= 6) "Likely pathogenic" else if (total >= 0) "Uncertain significance"
  else if (total >= -6) "Likely benign" else "Benign"
}

# All the criteria that can be chosen by hand
acmg_all_criteria <- c("PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2_Supporting", "PM2", "PM3", "PM4", "PM5", "PM6",
                       "PP1", "PP2", "PP3", "PP3_Moderate", "PP3_Strong", "PP4",
                       "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP4_Moderate", "BP4_Strong", "BP5", "BP7")

# Is this row de novo in the family (PED or trio with both parents)?
is_de_novo_row <- function(v, ped) {
  if (is.null(ped)) return(FALSE)
  row <- as.data.frame(v, stringsAsFactors = FALSE, check.names = FALSE)
  nrow(inheritance_filter(row, "De novo", ped)) == 1
}

# Column of AlphaMissense scores: added by PriorR or from the annotation of the file
alphamissense_column <- function(columns) {
  cols <- grep("alphamissense|am_pathogenicity", columns, ignore.case = TRUE, value = TRUE)
  cols <- cols[!grepl("class|pred$|rankscore", cols, ignore.case = TRUE)]
  if (length(cols)) cols[1] else NULL
}

# ---- Reanalysis ----
# What the annotation said when a variant was classified; compared later with the current annotation
snapshot_fields <- c("CLNSIG", "CLNREVSTAT", "gnomADg_AF_popmax", "gnomADe_AF_popmax", "Consequence")

variant_snapshot <- function(v) {
  out <- lapply(stats::setNames(snapshot_fields, snapshot_fields), function(f)
    if (f %in% names(v) && !is.na(v[[f]]) && nzchar(as.character(v[[f]]))) as.character(v[[f]]) else NULL)
  out <- Filter(Negate(is.null), out)
  # green PanelApp panels of the gene (in the PanelApp copy of the user)
  gene <- as.character(v$SYMBOL %||% NA)
  if (!is.na(gene) && nzchar(gene)) out$panelapp <- gene_panelapp_panels(gene)
  out$acmg_suggestion <- acmg_class(acmg_suggest(v)$criterion)
  out$date <- format(Sys.Date())
  out
}

# Ids of the PanelApp panels where the gene is green
gene_panelapp_panels <- function(gene, catalog = panel_catalog()) {
  pa <- catalog[startsWith(catalog$id, "panelapp_"), ]
  pa$id[vapply(pa$genes, function(g) gene %in% g, logical(1))]
}

# Classified variants of the table whose annotation changed since their classification:
# data.frame(variant_key, gene, classification, author, date, changes)
reanalysis <- function(df, comments) {
  empty <- data.frame(variant_key = character(0), gene = character(0), classification = character(0),
                      author = character(0), date = character(0), changes = character(0))
  if (!nrow(df) || !nrow(comments)) return(empty)
  classified <- comments[!is.na(comments$classification), ]
  classified <- classified[order(classified$created_at, decreasing = TRUE), ]
  classified <- classified[!duplicated(classified$variant_key), ]          # latest classification
  keys <- variant_key(df$CHROM, df$POS, df$REF, df$ALT)
  classified <- classified[classified$variant_key %in% keys, ]
  if (!nrow(classified)) return(empty)
  # one row per variant (the canonical transcript when there is one)
  o <- order(keys, if ("CANONICAL" %in% names(df)) df$CANONICAL != "YES" else 0)
  rows <- df[o, , drop = FALSE][!duplicated(keys[o]), , drop = FALSE]
  rows_keys <- variant_key(rows$CHROM, rows$POS, rows$REF, rows$ALT)
  freq_level <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!length(x) || is.na(x)) "absent" else if (x > 0.05) "> 0.05" else if (x > 0.01) "> 0.01"
    else if (x > 1e-4) "0.0001-0.01" else "<= 0.0001"
  }
  now <- function(v, f) if (f %in% names(v) && !is.na(v[[f]]) && nzchar(as.character(v[[f]]))) as.character(v[[f]]) else NA_character_
  changes <- lapply(seq_len(nrow(classified)), function(i) {
    c <- classified[i, ]
    v <- rows[match(c$variant_key, rows_keys), , drop = FALSE]
    snap <- c$snapshot[[1]]
    out <- character(0)
    if (!is.null(snap)) {
      before <- snap$CLNSIG %||% NA_character_
      if (!identical(before, now(v, "CLNSIG")))
        out <- c(out, sprintf("ClinVar: %s → %s", before %|NA|% "not in ClinVar", now(v, "CLNSIG") %|NA|% "not in ClinVar"))
      for (f in c("gnomADg_AF_popmax", "gnomADe_AF_popmax")) if (f %in% names(v)) {
        a <- freq_level(snap[[f]]); b <- freq_level(now(v, f))
        if (a != b) out <- c(out, sprintf("%s: %s → %s", f, a, b))
      }
      if (!is.null(snap$panelapp) && !is.na(c$gene)) {
        new <- setdiff(gene_panelapp_panels(c$gene), unlist(snap$panelapp))
        if (length(new)) {
          cat <- panel_catalog()
          out <- c(out, paste("New green PanelApp panels:", paste(cat$name[match(new, cat$id)], collapse = "; ")))
        }
      }
    }
    if (!is.null(snap$acmg_suggestion)) {
      suggested <- acmg_class(acmg_suggest(v)$criterion)
      if (suggested != snap$acmg_suggestion) out <- c(out, sprintf("ACMG suggestion: %s → %s", snap$acmg_suggestion, suggested))
    }
    paste(out, collapse = " | ")
  })
  out <- data.frame(variant_key = classified$variant_key, gene = classified$gene, classification = classified$classification,
                    author = classified$author, date = substr(classified$created_at, 1, 10), changes = unlist(changes),
                    stringsAsFactors = FALSE)
  out[nzchar(out$changes), , drop = FALSE]
}

`%|NA|%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
