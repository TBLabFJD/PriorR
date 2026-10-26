# Details of the selected variant: genotypes, frequencies, predictors, OMIM and comments.

variant_info_ui <- function(id) {
  ns <- NS(id)
  shinydashboardPlus::box(
    title = "Variant information", width = NULL, status = "success", collapsible = TRUE, collapsed = FALSE,
    info_css,
    uiOutput(ns("empty")),
    div(id = ns("details"),
      uiOutput(ns("header")),
      fluidRow(
        column(4,
               card("Transcript", "dna", "#00897b", uiOutput(ns("summary"))),
               card("Links", "external-link-alt", "#5e35b1", uiOutput(ns("links")))),
        column(4,
               card("Genotypes", "users", "#1565c0", uiOutput(ns("genotypes"))),
               card("Read depth (DP)", "chart-bar", "#0277bd", uiOutput(ns("depth_text")), plotOutput(ns("depth"), height = 100))),
        column(4,
               card("Population frequency", "globe-europe", "#2e7d32", uiOutput(ns("frequencies"))))),
      # full width: predictors on the left, scores in two columns on the right
      card("Pathogenicity", "exclamation-triangle", "#e65100",
           fluidRow(column(4, uiOutput(ns("predictors"))), column(8, uiOutput(ns("scores"))))),
      fluidRow(
        column(4, card("OMIM phenotype", "notes-medical", "#ad1457", uiOutput(ns("omim")))),
        column(4, card("ACMG (suggested)", "balance-scale", "#283593", uiOutput(ns("acmg")))),
        column(4, card("Comments", "comments", "#6d4c41", comments_ui(ns("comments")))))))
}

variant_label <- function(v) {
  allele <- function(x) if (is.na(x) || x == "") "-" else as.character(x)
  paste(allele(v$SYMBOL), v$CHROM, v$POS, allele(v$REF), ">", allele(v$ALT))
}

# variant: reactive with the selected row (all columns) or NULL; data: all the loaded variants;
# ped: the family (PED or trio) for the de novo evidence of ACMG
variant_info_server <- function(id, variant, data, user, ped = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    value <- function(col) {
      v <- variant()
      if (is.null(v) || !col %in% names(v) || is.na(v[[col]])) "-" else as.character(v[[col]])
    }
    lines <- function(...) {
      items <- list(...)
      HTML(paste(sprintf("<b>%s</b> %s", names(items), vapply(items, value, "")), collapse = "<br/>"))
    }

    output$empty <- renderUI(if (is.null(variant())) p(class = "text-muted", "Select a variant in the table."))
    observe(shinyjs::toggle("details", condition = !is.null(variant())))
    output$header <- renderUI({
      v <- req(variant())
      # IMPACT of VEP when the file has it; otherwise from the consequence (style.R)
      impact <- if (value("IMPACT") %in% names(impact_color)) value("IMPACT") else consequence_impact(value("Consequence"))
      priority <- priority_class(as.data.frame(v))
      div(class = "pr-header",
          span(class = "pr-gene", value("SYMBOL")),
          span(class = "pr-pos", sprintf("%s:%s %s>%s", value("CHROM"), value("POS"), value("REF"), value("ALT"))),
          div(style = "margin-top: 6px;",
              badge(paste(impact, "impact"), impact_color[[impact]]),
              badge(gsub("_", " ", strsplit(value("Consequence"), ",")[[1]][1]), impact_color[[impact]]),
              if (value("CLNSIG") != "-") badge(paste("ClinVar:", gsub("_", " ", value("CLNSIG"))), clinvar_color(value("CLNSIG"))),
              if (!is.na(priority)) badge(paste("PriorR", priority), if (priority == "red") "#c62828" else "#ef6c00"),
              if (value("VARIANT_CLASS") != "-") badge(value("VARIANT_CLASS"), "#546e7a"),
              if (value("CANONICAL") == "YES") badge("canonical", "#00897b")))
    })
    # External links (Program_files/links.tsv)
    output$links <- renderUI({
      v <- req(variant())
      as_buttons <- function(l) lapply(seq_len(nrow(l)), function(i)
        tags$a(href = l$url[i], target = "_blank", title = l$description[i], class = "btn btn-default btn-xs",
               style = "margin: 0 4px 4px 0;", l$name[i]))
      gene <- variant_links(v, "gene")
      variant <- variant_links(v, "variant")
      as_buttons <- function(l, class) lapply(seq_len(nrow(l)), function(i)
        tags$a(href = l$url[i], target = "_blank", title = l$description[i], class = paste("btn btn-xs", class), l$name[i]))
      div(class = "pr-links",
        if (nrow(gene)) div(strong(paste0(v$SYMBOL, " ")), as_buttons(gene, "pr-link-gene")),
        if (nrow(variant)) div(style = "margin-top: 4px;", strong("Variant "), as_buttons(variant, "pr-link-variant")))
    })

    output$summary <- renderUI({
      req(variant())
      tagList(kv("Consequence", gsub("_", " ", value("Consequence"))), kv("HGVSc", value("HGVSc")),
              kv("HGVSp", value("HGVSp")), kv("Exon / intron", paste(value("EXON"), "/", value("INTRON"))),
              kv("Feature", value("Feature")), kv("Existing variation", value("Existing_variation")))
    })

    # One row per sample: genotype, depth, alleles, VAF and the genotype of each variant caller
    output$genotypes <- renderUI({
      v <- req(variant())
      samples <- sample_names(names(v))
      req(length(samples) > 0)
      cell <- function(col) if (col %in% names(v) && !is.na(v[[col]])) as.character(v[[col]]) else "-"
      tags$table(class = "table table-condensed", style = "font-size: 13px; margin-bottom: 0;",
        tags$thead(tags$tr(lapply(c("Sample", "GT", "DP", "AD", "VAF", "Callers"), tags$th))),
        tags$tbody(lapply(samples, function(s) tags$tr(
          tags$td(strong(s)), tags$td(genotype_chip(cell(paste0(s, "_GT")))), tags$td(cell(paste0(s, "_DP"))),
          tags$td(cell(paste0(s, "_AD"))), tags$td(cell(paste0(s, "_VAF"))),
          tags$td(lapply(sample_callers(names(v), s), function(c)
            span(style = "margin-right: 4px;", span(class = "text-muted", c), genotype_chip(cell(paste0(s, "_", c, "_GT"))))))))))
    })
    # Frequencies: green rare, orange low frequency, grey common
    output$frequencies <- renderUI({
      req(variant())
      freq <- function(label, col) {
        x <- suppressWarnings(as.numeric(value(col)))
        colors <- if (is.na(x)) c(bg = "#e8f5e9", fg = "#2e7d32") else if (x < 0.001) c(bg = "#e8f5e9", fg = "#2e7d32")
                  else if (x < 0.01) c(bg = "#fff3e0", fg = "#e65100") else c(bg = "#eceff1", fg = "#546e7a")
        kv(label, chip(if (is.na(x)) "absent" else format(signif(x, 3), scientific = x < 1e-4), colors))
      }
      tagList(freq("gnomAD genomes", "gnomADg_AF"), freq(paste0("gnomAD popmax (", value("gnomADg_popmax"), ")"), "gnomADg_AF_popmax"),
              freq("gnomAD exomes", "gnomADe_AF"), freq("FJD", "FJD_MAF_AF"), freq("Spanish (CSVS)", "CSVS_AF"),
              helpText(chip("< 0.001", c(bg = "#e8f5e9", fg = "#2e7d32")), chip("< 0.01", c(bg = "#fff3e0", fg = "#e65100")),
                       chip("common", c(bg = "#eceff1", fg = "#546e7a"))))
    })

    # Scores as colored bars: the longer and redder, the more pathogenic it looks
    output$scores <- renderUI({
      req(variant())
      gauge <- function(label, col, max, thresholds) {
        x <- suppressWarnings(as.numeric(value(col)))
        colors <- level_color(x, thresholds)
        width <- if (is.na(x)) 0 else min(100, max(3, 100 * x / max))
        div(style = "display: flex; align-items: center; margin: 3px 0; font-size: 13px;",
            span(style = "width: 95px; color: #607d8b;", label),
            div(style = "flex: 1; height: 10px; background: #eceff1; border-radius: 5px; margin: 0 8px; position: relative;",
                div(style = sprintf("width: %s%%; height: 10px; border-radius: 5px; background: %s;", width, colors[["fg"]])),
                lapply(thresholds, function(t) div(style = sprintf("position: absolute; left: %s%%; top: -2px; height: 14px; border-left: 1px dashed #b0bec5;", min(100, 100 * t / max))))),
            span(style = sprintf("width: 52px; text-align: right; font-weight: 600; color: %s;", colors[["fg"]]),
                 if (is.na(x)) "-" else format(round(x, 3), nsmall = 0)))
      }
      am <- alphamissense_column(names(variant()))
      tagList(div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); column-gap: 28px;",
                  gauge("CADD PHRED", "CADD_PHRED", 50, c(15, 20, 25)), gauge("REVEL", "REVELScore", 1, c(0.3, 0.5, 0.75)),
                  if (!is.null(am)) gauge("AlphaMissense", am, 1, c(0.34, 0.564, 0.8)),
                  gauge("SpliceAI max", "SpliceAI_DS_Max", 1, c(0.1, 0.2, 0.5)), gauge("MutScore", "MutScore", 1, c(0.5, 0.65, 0.8)),
                  gauge("ADA score", "ada_score", 1, c(0.4, 0.6, 0.8)), gauge("RF score", "rf_score", 1, c(0.4, 0.6, 0.8))),
              if (value("gnomAD_LOEUF") != "-")
                div(style = "margin-top: 6px; font-size: 13px;", span(style = "color: #607d8b;", "Gene constraint (gnomAD v4.1): "),
                    chip(paste("LOEUF", value("gnomAD_LOEUF")), level_color(value("gnomAD_LOEUF"), c(0.6, 0.35, 0.2), high_is_bad = FALSE)),
                    chip(paste("pLI", value("gnomAD_pLI")), level_color(value("gnomAD_pLI"), c(0.5, 0.9, 0.99))),
                    chip(paste("mis Z", value("gnomAD_mis_z")), level_color(value("gnomAD_mis_z"), c(2, 3.09, 4)))),
              helpText("Dashed lines: thresholds. Red: likely damaging. LOEUF < 0.6: intolerant to loss of function."))
    })

    # ACMG criteria suggested by the annotation (acmg.R)
    acmg <- reactive({
      v <- req(variant())
      acmg_suggest(v, de_novo = is_de_novo_row(v, ped()))
    })
    output$acmg <- renderUI({
      a <- acmg()
      cls <- acmg_class(a$criterion)
      cls_color <- clinvar_color(cls)
      tagList(
        div(style = "margin-bottom: 6px;", badge(cls, cls_color), span(class = "text-muted", sprintf(" %+d points", acmg_total(a$criterion)))),
        if (!nrow(a)) p(class = "text-muted", "No criteria from the annotation.")
        else tags$ul(class = "list-unstyled", style = "font-size: 13px;", lapply(seq_len(nrow(a)), function(i)
          tags$li(style = "margin-bottom: 3px;", badge(a$criterion[i], if (a$points[i] > 0) "#c62828" else "#2e7d32"), a$reason[i]))),
        helpText("Only from the annotation (frequency, consequence, REVEL/CADD, SpliceAI, de novo): add the rest (segregation, functional, phenotype...) when classifying."))
    })

    output$omim <- renderUI({
      req(variant())
      phenotypes <- strsplit(value("OMIM_phenotype"), ";", fixed = TRUE)[[1]]
      tagList(lapply(trimws(phenotypes), function(ph) div(style = "padding: 2px 0;", inheritance_badges(ph), ph)))
    })

    # Depth of this variant compared with all the variants of the file, per sample
    depth_data <- reactive({
      v <- req(variant())
      df <- req(data())
      samples <- Filter(function(s) paste0(s, "_DP") %in% names(df), sample_names(names(df)))
      req(length(samples) > 0)
      variants <- df[!duplicated(df[, c("CHROM", "POS")]), , drop = FALSE]
      lapply(stats::setNames(samples, samples), function(s) {
        all <- suppressWarnings(as.numeric(variants[[paste0(s, "_DP")]]))
        all <- all[!is.na(all)]
        this <- suppressWarnings(as.numeric(v[[paste0(s, "_DP")]]))
        list(all = all, this = this, percentile = if (is.na(this) || !length(all)) NA else round(100 * mean(all <= this)))
      })
    })

    output$depth_text <- renderUI({
      d <- depth_data()
      tagList(lapply(names(d), function(s) {
        x <- d[[s]]
        if (is.na(x$this)) return(p(s, ": no depth"))
        verdict <- if (x$this < 10 || x$percentile <= 5) span(class = "text-danger", strong("low"))
                   else if (x$percentile >= 95) span(class = "text-warning", strong("high"))
                   else span(class = "text-success", "normal")
        p(style = "margin-bottom: 2px;", sprintf("%s: DP %s, more than %s%% of the variants of the file \u2014 ", s, x$this, x$percentile), verdict)
      }), helpText("Box: depth of all the variants of the file. Red dot: this variant."))
    })

    output$depth <- renderPlot({
      d <- depth_data()
      all <- do.call(rbind, lapply(names(d), function(s) data.frame(sample = s, DP = d[[s]]$all)))
      this <- do.call(rbind, lapply(names(d), function(s) data.frame(sample = s, DP = d[[s]]$this)))
      limit <- max(stats::quantile(all$DP, 0.99, na.rm = TRUE), this$DP, na.rm = TRUE)
      ggplot2::ggplot(all, ggplot2::aes(x = DP, y = sample)) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#e8f5e9", color = "#5f7f6a", width = 0.5) +
        ggplot2::geom_point(data = this, color = "red", size = 3) +
        ggplot2::coord_cartesian(xlim = c(0, limit)) +
        ggplot2::labs(x = "Read depth (DP)", y = NULL) +
        ggplot2::theme_minimal(base_size = 12)
    })

    # Predictors (dbNSFP, only for missense variants): a bar and the name of each predictor
    output$predictors <- renderUI({
      v <- req(variant())
      counts <- suppressWarnings(as.numeric(c(v$N_Pathogenic_pred, v$N_Benign_pred)))
      if (any(is.na(counts)) || sum(counts) == 0) {
        consequence <- gsub("_", " ", strsplit(value("Consequence"), ",")[[1]][1])
        message <- if (grepl("missense|protein altering|inframe", consequence))
          sprintf("No predictor results in the annotation for this %s (dbNSFP has no prediction for this transcript or variant). See the scores below.", consequence)
        else sprintf("No missense predictors: they only apply to amino acid changes (this one is %s). See the scores below.", consequence)
        return(p(class = "text-muted", style = "font-size: 13px;", icon("info-circle"), message))
      }
      pct <- round(100 * counts / sum(counts))
      names_of <- function(col) {
        x <- value(col)
        if (x == "-") character(0) else gsub("_pred$", "", trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
      }
      tagList(
        div(style = "display: flex; height: 22px; border-radius: 11px; overflow: hidden; margin: 4px 0 2px 0; font-size: 12px; color: white; font-weight: 600;",
            div(style = sprintf("width: %s%%; background: #e64a19; text-align: center; line-height: 22px;", pct[1]), if (pct[1] >= 15) counts[1]),
            div(style = sprintf("width: %s%%; background: #43a047; text-align: center; line-height: 22px;", pct[2]), if (pct[2] >= 15) counts[2])),
        div(style = "font-size: 12px; margin-bottom: 4px;", span(style = "color: #e64a19; font-weight: 600;", counts[1], "pathogenic"),
            " \u00b7 ", span(style = "color: #43a047; font-weight: 600;", counts[2], "benign"), " predictors"),
        div(style = "margin-bottom: 8px; line-height: 1.9;",
            lapply(names_of("Pathogenic_pred"), function(n) chip(n, c(bg = "#fbe9e7", fg = "#d84315"))),
            lapply(names_of("Benign_pred"), function(n) chip(n, c(bg = "#e8f5e9", fg = "#2e7d32")))))
    })

    selected <- reactive({
      v <- variant()
      if (is.null(v)) return(NULL)
      list(key = variant_key(v$CHROM, v$POS, v$REF, v$ALT), gene = as.character(v$SYMBOL), label = variant_label(v),
           acmg = acmg()$criterion, snapshot = variant_snapshot(v))
    })
    comments_server("comments", user, selected)
  })
}
