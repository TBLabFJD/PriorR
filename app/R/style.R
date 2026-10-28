# Visual helpers: colored badges and cards of the variant details.

info_css <- shiny::tags$style(shiny::HTML("
  .pr-header { background: linear-gradient(90deg, #e8f5e9, #ffffff); border-radius: 6px; padding: 10px 14px; margin-bottom: 12px; }
  .pr-header .pr-gene { font-size: 26px; font-weight: 700; color: #1b5e20; margin-right: 10px; }
  .pr-header .pr-pos { font-size: 16px; color: #37474f; font-family: monospace; }
  .pr-badge { display: inline-block; padding: 3px 9px; border-radius: 12px; font-size: 12px; font-weight: 600;
              margin: 2px 4px 2px 0; color: #fff; white-space: nowrap; }
  .pr-card { background: #fff; border: 1px solid #e3e7ea; border-left: 4px solid var(--card-color, #00a65a);
             border-radius: 6px; padding: 8px 12px; margin-bottom: 12px; }
  .pr-card h5 { margin-top: 2px; font-weight: 700; color: var(--card-color, #00a65a); }
  .pr-kv { display: flex; justify-content: space-between; padding: 2px 0; border-bottom: 1px dotted #eceff1; }
  .pr-kv:last-child { border-bottom: none; }
  .pr-kv .pr-k { color: #607d8b; }
  .pr-chip { display: inline-block; min-width: 54px; text-align: center; padding: 1px 8px; border-radius: 10px;
             font-weight: 600; font-size: 12px; }
  .pr-links a { margin: 0 4px 4px 0; }
  .pr-link-gene { background: #e0f2f1; border-color: #80cbc4; color: #00695c; }
  .pr-link-variant { background: #ede7f6; border-color: #b39ddb; color: #4527a0; }
"))

badge <- function(text, color) span(class = "pr-badge", style = sprintf("background: %s;", color), text)

# Color of a value between thresholds: c(low, mid, high) -> grey / yellow / orange / red
level_color <- function(x, thresholds, high_is_bad = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  if (is.na(x)) return(c(bg = "#eceff1", fg = "#90a4ae"))
  level <- sum(if (high_is_bad) x >= thresholds else x <= thresholds)
  palette <- list(c(bg = "#eceff1", fg = "#546e7a"), c(bg = "#fff8e1", fg = "#f57f17"),
                  c(bg = "#ffe0b2", fg = "#e65100"), c(bg = "#ffcdd2", fg = "#b71c1c"))
  palette[[level + 1]]
}

chip <- function(text, colors) span(class = "pr-chip", style = sprintf("background: %s; color: %s;", colors[["bg"]], colors[["fg"]]), text)

card <- function(title, icon_name, color, ...) {
  div(class = "pr-card", style = sprintf("--card-color: %s;", color), h5(icon(icon_name), title), ...)
}

kv <- function(key, value) div(class = "pr-kv", span(class = "pr-k", key), span(value))

# Impact of VEP consequences (vectorized)
consequence_impact <- function(consequence) {
  c <- tolower(as.character(consequence %||% ""))
  c[is.na(c)] <- ""
  ifelse(grepl("stop_gained|frameshift|splice_acceptor|splice_donor_variant|start_lost|stop_lost|transcript_ablation", c), "HIGH",
  ifelse(grepl("missense|inframe|protein_altering|splice_donor_5th|splice_donor_region", c), "MODERATE",
  ifelse(grepl("synonymous|splice_region|splice_polypyrimidine|stop_retained|start_retained", c), "LOW", "MODIFIER")))
}
impact_color <- c(HIGH = "#d32f2f", MODERATE = "#ef6c00", LOW = "#f9a825", MODIFIER = "#78909c")

clinvar_color <- function(clnsig) {
  s <- tolower(clnsig %||% "")
  if (grepl("conflicting", s)) "#7b1fa2"
  else if (grepl("likely_pathogenic|likely pathogenic", s) && !grepl("^pathogenic", s)) "#ef6c00"
  else if (grepl("pathogenic", s)) "#c62828"
  else if (grepl("uncertain", s)) "#f9a825"
  else if (grepl("benign", s)) "#2e7d32"
  else "#78909c"
}

genotype_chip <- function(gt) {
  colors <- if (gt_hom_alt(gt)) c(bg = "#ede7f6", fg = "#5e35b1") else if (gt_het(gt)) c(bg = "#e3f2fd", fg = "#1565c0")
            else if (gt_ref(gt)) c(bg = "#e8f5e9", fg = "#2e7d32") else c(bg = "#eceff1", fg = "#90a4ae")
  label <- if (gt_hom_alt(gt)) "hom" else if (gt_het(gt)) "het" else if (gt_ref(gt)) "ref" else "no call"
  chip(sprintf("%s %s", gt, label), colors)
}

inheritance_badges <- function(omim) {
  o <- tolower(omim %||% "")
  tagList(if (grepl("autosomal recessive", o)) badge("AR", "#5e35b1"),
          if (grepl("autosomal dominant", o)) badge("AD", "#00838f"),
          if (grepl("x-linked", o)) badge("XL", "#ad1457"),
          if (grepl("mitochondrial", o)) badge("MT", "#6d4c41"))
}
