# Boxes above the SNV table that explain the PriorR priority and GLOWgenes ranking views.

snv_view_info_server <- function(input, output, session, data, filtered, glow_ready, glow_panels) {
  ns <- session$ns
  # PriorR priority view: how many red / orange variants and what they mean (prioritization.R)
  output$priority_info <- renderUI({
    req(identical(input$view, "priority"), data())
    DF <- filtered()
    cls <- attr(DF, "priority") %||% priority_class(DF)
    keys <- variant_key(DF$CHROM, DF$POS, DF$REF, DF$ALT)
    n <- function(c) length(unique(keys[cls %in% c]))
    rules <- function(x, color) tags$ul(style = "margin: 2px 0 6px 0;", lapply(x, function(r) tags$li(r)))
    div(class = "alert", style = "padding: 8px 12px; margin-bottom: 8px; background: #fff3e0; border: 1px solid #ffcc80; color: #4e342e;",
        strong(sprintf("PriorR priority: %s red and %s orange variants first, then the rest (with the sidebar filters). ",
                       format(n("red"), big.mark = ","), format(n("orange"), big.mark = ","))),
        tags$details(tags$summary(style = "cursor: pointer;", "What red and orange mean"),
          span(style = "color: #c62828; font-weight: 600;", "Red"), paste0(" (any of them; ", priority_rules$red_never, "):"),
          rules(priority_rules$red),
          span(style = "color: #ef6c00; font-weight: 600;", "Orange"), " (any of them):",
          rules(priority_rules$orange)))
  })

  output$glow_info <- renderUI({
    req(identical(input$view, "glow"), data())
    if (!glow_ready()) return(div(class = "alert alert-warning", style = "padding: 8px;",
      if (!has_glow_column(data())) "This file has no GLOWgenes column: choose a precomputed ranking or run a new one in the sidebar (GLOWgenes ranking)."
      else "Choose a precomputed ranking or run a new one in the sidebar (GLOWgenes ranking)."))
    out <- filtered()
    s <- glow_settings()
    e <- s$evidence
    source <- if (is.null(glow_panels())) "the GLOWgenes column of the file" else paste("the rankings", paste(gsub("_", " ", sub("^[^:]*:", "", glow_panels())), collapse = ", "))
    evidence <- c(if (isTRUE(e$clinvar_pathogenic)) "ClinVar pathogenic", if (isTRUE(e$high_impact)) "high impact",
                  if (!is.null(e$min_cadd)) paste("CADD \u2265", e$min_cadd), if (!is.null(e$min_revel)) paste("REVEL \u2265", e$min_revel),
                  if (!is.null(e$min_spliceai)) paste("SpliceAI \u2265", e$min_spliceai), if (isTRUE(e$priorr_priority)) "PriorR red / orange")
    div(class = "alert alert-success", style = "padding: 8px; margin-bottom: 8px;",
        strong(sprintf("GLOWgenes ranking: %s variants in %s genes, ordered by %s. ", format(nrow(out), big.mark = ","),
                       length(unique(stats::na.omit(out$SYMBOL))), source)),
        sprintf("Rare (frequency \u2264 %s)%s and at least one of: %s. Rank 0 = genes of the disease panel%s.",
                max(unlist(s$max_af)), if (isTRUE(s$exclude_clinvar_benign)) ", not ClinVar benign" else "", paste(evidence, collapse = ", "),
                if (isTRUE(input$glow_hide_rank0)) " (hidden: only new candidates)" else ""))
  })
}
