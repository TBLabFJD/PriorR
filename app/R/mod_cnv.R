# CNV analysis (AnnotSV files): sidebar with the filters, table, details and download.
# By default nothing is filtered (one row per CNV: Annotation mode "full").

cnv_sidebar_ui <- function(id) {
  ns <- NS(id)
  picker <- function(id, label, choices = NULL, selected = choices)
    shinyWidgets::pickerInput(ns(id), label, choices = choices, selected = selected,
                              options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple = TRUE)
  check <- function(id, label, value = FALSE)
    shinyWidgets::prettyCheckbox(ns(id), label, value = value, outline = TRUE, status = "success")
  panels <- function(id, label, specific = NA) selectizeInput(ns(id), label, choices = panel_choices(specific = specific), multiple = TRUE,
                                                              options = list(placeholder = "Search a panel", maxOptions = 2000))
  shinydashboard::sidebarMenu(
    id = ns("menu"),
    shinydashboard::menuItem("Load CNVs", icon = icon("upload"), startExpanded = TRUE,
      open_file_ui(ns("file"), "Open CNV file", "Open a CNV file of AnnotSV (.tsv)"),
      picker("samples", "Samples")),
    shinydashboard::menuItem("Columns", icon = icon("columns"),
      picker("columns", NULL)),
    shinydashboard::menuItem("Filters", icon = icon("ruler-horizontal"),
      picker("annotation_mode", "Annotation mode (full: one row per CNV; split: one per gene)", c("full", "split"), "full"),
      picker("type", "Type"),
      picker("region", "Location"),
      shinyWidgets::sliderTextInput(ns("min_acmg"), "Min ACMG class", choices = 1:5, selected = 1),
      div(id = ns("programs_box"), sliderInput(ns("min_programs"), "Min number of programs", min = 1, max = 3, value = 1, step = 1)),
      check("order_by_class", "Order by ACMG class and AnnotSV score", TRUE)),
    shinydashboard::menuItem("Genes and panels", icon = icon("dna"),
      textInput(ns("genes"), "Your genes (comma separated)"),
      panels("panels", "Gene panels (PanelApp, RetiGene, Orphanet, saved)", specific = FALSE),
      check("amber", "Include amber genes of PanelApp / RetiGene candidates"),
      check("red", "Include amber and red genes of PanelApp"),
      panels("specific_panels", "Disease-specific genes (Benítez et al., 2026)", specific = TRUE),
      fileInput(ns("panel_file"), "...or upload a panel (one gene per line)", accept = c(".txt", ".csv", ".tsv")),
      panels("excluded_panels", "Exclude the genes of these panels")),
    shinydashboard::menuItem("GLOWgenes ranking", icon = icon("sort-amount-down"),
      selectizeInput(ns("glow_panels"), "Rankings (the best rank of the genes of each CNV)", choices = glow_rank_choices(), multiple = TRUE,
                     options = list(placeholder = "Orphanet, PanelApp or computed by the lab", maxOptions = 1000)),
      check("glow_order", "Order by GLOWgenes rank"),
      check("glow_limit", "Only CNVs with a gene up to a GLOWgenes rank"),
      conditionalPanel("input.glow_limit", ns = ns,
        sliderInput(ns("glow"), "Max GLOWgenes rank", min = 0, max = 1800, value = 100))),
    shinydashboard::menuItem("Phenotypes / HPO", icon = icon("user-cog"),
      selectizeInput(ns("phenotype"), "HPO terms", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Type a phenotype or HP: identifier")),
      textInput(ns("hpo_ids"), "...or paste HPO identifiers", placeholder = "HP:0000510, HP:0000556"),
      helpText("The genes of the terms (and of their subterms) work like a gene panel.")),
    shinydashboard::menuItem("Frequencies", icon = icon("users"),
      sliderInput(ns("max_gain_af"), "Max gain AF (B_gain_AFmax)", min = 0, max = 1, value = 1),
      sliderInput(ns("max_loss_af"), "Max loss AF (B_loss_AFmax)", min = 0, max = 1, value = 1)),
    shinydashboard::menuItem("Save data", icon = icon("save"),
      downloadButton(ns("download"), "Download filtered data", class = "btn-primary", style = "margin: 6px 15px;")))
}

cnv_main_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("loaded_info")),
    DT::DTOutput(ns("table")),
    br(),
    shinydashboardPlus::box(
      title = "CNV information", width = NULL, status = "success", collapsible = TRUE,
      uiOutput(ns("details"))),
    reads_ui(ns("reads")))
}

cnv_default_columns <- function(df) {
  wanted <- c("GLOW_rank", "ACMG_class", "AnnotSV_ranking_score", "Samples_ID", paste0(cnv_samples(df), "_GT"),
              "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "N_PROGRAMS", "Annotation_mode",
              "CytoBand", "Gene_name", "Gene_count", "Tx", "Exon_count", "Location", "Location2", "Frameshift",
              "P_gain_phen", "P_loss_phen", "B_gain_AFmax", "B_loss_AFmax", "GenCC_disease", "OMIM_phenotype",
              "OMIM_inheritance", "LOEUF_bin", "GnomAD_pLI", "DDD_disease")
  intersect(wanted, names(df))
}

cnv_server <- function(id, user) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    cnv_file <- open_file_server("file", input, output)
    # the PanelApp copy of the user (panels.R)
    observeEvent(user(), {
      updateSelectizeInput(session, "panels", choices = panel_choices(specific = FALSE), selected = isolate(input$panels))
      updateSelectizeInput(session, "excluded_panels", choices = panel_choices(), selected = isolate(input$excluded_panels))
    })
    updateSelectizeInput(session, "phenotype", choices = hpo_choices(), server = TRUE)
    observeEvent(input$hpo_ids, {
      ids <- parse_hpo_ids(input$hpo_ids)
      req(length(ids) > 0)
      known <- ids[ids %in% hpo_choices()]
      if (length(known)) updateSelectizeInput(session, "phenotype", choices = hpo_choices(),
                                              selected = unique(c(input$phenotype, known)), server = TRUE)
      if (length(known) < length(ids))
        showNotification(paste("HPO terms without genes or unknown:", paste(setdiff(ids, known), collapse = ", ")), type = "warning")
      updateTextInput(session, "hpo_ids", value = "")
    }, ignoreInit = TRUE)
    data <- reactive({
      req(cnv_file(), user())
      df <- withProgress(message = "Loading CNVs...", value = 0.5, tryCatch(read_cnvs(cnv_file()), error = function(e) NULL))
      validate(need(!is.null(df) && all(c("SV_chrom", "SV_start", "SV_end") %in% names(df)),
                    paste("This file can't be read as a CNV table: PriorR opens the tab-separated files of AnnotSV (SV_chrom, SV_start, SV_end, ...):",
                          basename(cnv_file()))))
      df
    })

    observeEvent(data(), {
      df <- data()
      choices_of <- function(col) if (col %in% names(df)) sort(unique(as.character(df[[col]]))) else character(0)
      samples <- sort(unique(trimws(unlist(strsplit(stats::na.omit(as.character(df$Samples_ID)), ",", fixed = TRUE)))))
      shinyWidgets::updatePickerInput(session, "samples", choices = samples, selected = samples)
      shinyWidgets::updatePickerInput(session, "columns", choices = column_choices(names(df)), selected = cnv_default_columns(df))
      shinyWidgets::updatePickerInput(session, "region", choices = c(choices_of("Location2"), "NA"), selected = c(choices_of("Location2"), "NA"))
      shinyWidgets::updatePickerInput(session, "type", choices = choices_of("SV_type"), selected = choices_of("SV_type"))
      # "Number of programs" only for files that have it
      shinyjs::toggle("programs_box", condition = "N_PROGRAMS" %in% names(df))
    })

    filtered <- reactive({
      df <- data()
      panel <- c(panel_genes(c(input$panels, input$specific_panels), amber = isTRUE(input$amber), red = isTRUE(input$red)),
                 if (!is.null(input$panel_file)) read_gene_list(input$panel_file$datapath), hpo_genes(input$phenotype))
      excluded <- panel_genes(input$excluded_panels, amber = TRUE)
      all_samples <- length(input$samples) == 0 || setequal(input$samples, unique(trimws(unlist(strsplit(stats::na.omit(as.character(df$Samples_ID)), ",")))))
      filter_cnv(df, cnv_params(
        samples = if (!all_samples) input$samples, annotation_mode = input$annotation_mode, type = input$type,
        region = input$region, min_acmg = input$min_acmg, min_programs = input$min_programs,
        genes = split_genes(input$genes), panel_genes = if (length(panel)) unique(panel),
        excluded_genes = if (length(excluded)) excluded, max_gain_af = input$max_gain_af, max_loss_af = input$max_loss_af,
        order_by_class = input$order_by_class, columns = input$columns,
        glow_sources = input$glow_panels, glow_max = if (isTRUE(input$glow_limit) && length(input$glow_panels)) input$glow,
        glow_order = isTRUE(input$glow_order)))
    })

    output$loaded_info <- renderUI({
      req(data())
      p(class = "text-muted", sprintf("%s rows loaded · %s shown", format(nrow(data()), big.mark = ","),
                                      format(nrow(filtered()), big.mark = ",")))
    })

    # Displayed columns: the key ones first (ACMG class, score, genotype, position...)
    display_columns <- reactive({
      DF <- filtered()
      cols <- attr(DF, "display_columns")
      key <- cnv_default_columns(DF)
      c(intersect(key, cols), setdiff(cols, key))
    })

    output$table <- DT::renderDT({
      DF <- filtered()
      shown <- DF[, display_columns(), drop = FALSE]
      dt <- priorr_datatable(shown, ns("full_cell"))
      if ("ACMG_class" %in% names(shown))
        dt <- DT::formatStyle(dt, "ACMG_class", backgroundColor = DT::styleEqual(
          c("1", "2", "3", "4", "5"), c("#c8e6c9", "#e6f4c8", "#fff3b0", "#fdcf94", "#f7a1a1")))
      dt
    }, server = TRUE)
    observeEvent(input$full_cell, show_full_cell(input$full_cell))

    selected <- reactive({
      i <- input$table_rows_selected
      if (length(i) != 1 || i > nrow(filtered())) NULL else filtered()[i, , drop = FALSE]
    })
    # Reads of the selected CNV (the whole CNV or its breakpoints)
    reads_server("reads", selected, data, locus_of = cnv_locus)

    # Details of the selected CNV
    output$details <- renderUI({
      cnv <- selected()
      if (is.null(cnv)) return(p(class = "text-muted", "Select a CNV in the table."))
      value <- function(col) if (col %in% names(cnv) && !is.na(cnv[[col]]) && nzchar(as.character(cnv[[col]]))) as.character(cnv[[col]]) else "-"
      size <- abs(suppressWarnings(as.numeric(cnv$SV_length)))
      size_text <- if (is.na(size)) "-" else if (size >= 1e6) sprintf("%.1f Mb", size / 1e6) else if (size >= 1e3) sprintf("%.1f kb", size / 1e3) else sprintf("%d bp", as.integer(size))
      genes <- strsplit(value("Gene_name"), ";", fixed = TRUE)[[1]]
      genes <- genes[genes != "-"]
      samples <- cnv_samples(cnv)
      links <- variant_links(cnv, "cnv")
      acmg <- value("ACMG_class")
      acmg_color <- c("1" = "#2e7d32", "2" = "#7cb342", "3" = "#f9a825", "4" = "#ef6c00", "5" = "#c62828")[acmg]
      acmg_name <- c("1" = "benign", "2" = "likely benign", "3" = "uncertain", "4" = "likely pathogenic", "5" = "pathogenic")[acmg]
      type_color <- c(DEL = "#c62828", DUP = "#1565c0", INS = "#00838f", INV = "#6a1b9a", TRA = "#4e342e")[value("SV_type")]
      long_text <- function(col) div(style = "max-height: 110px; overflow-y: auto; word-break: break-word; font-size: 13px;",
                                     if (value(col) == "-") span(class = "text-muted", "-")
                                     else lapply(strsplit(value(col), ";", fixed = TRUE)[[1]], function(x) div(inheritance_badges(x), x)))
      tagList(
        info_css,
        div(class = "pr-header",
            span(class = "pr-gene", sprintf("chr%s:%s-%s", sub("^chr", "", value("SV_chrom")),
                                            format(as.numeric(cnv$SV_start), big.mark = ","), format(as.numeric(cnv$SV_end), big.mark = ","))),
            span(class = "pr-pos", value("CytoBand")),
            div(style = "margin-top: 6px;",
                badge(value("SV_type"), if (is.na(type_color)) "#546e7a" else type_color),
                badge(size_text, "#455a64"),
                if (!is.na(acmg_color)) badge(sprintf("ACMG class %s: %s", acmg, acmg_name), acmg_color),
                if (value("AnnotSV_ranking_score") != "-") badge(paste("AnnotSV score", value("AnnotSV_ranking_score")), "#546e7a"),
                badge(sprintf("%d genes", length(genes)), "#00897b"),
                lapply(samples, function(s) badge(paste(s, value(paste0(s, "_GT"))), "#1565c0")))),
        fluidRow(
          column(4,
                 card("CNV", "ruler-horizontal", "#00897b",
                      kv("Location", value("Location")), kv("Location (UTR / CDS)", value("Location2")),
                      kv("Exons", value("Exon_count")), kv("Frameshift", value("Frameshift")),
                      kv("Benign gain AF max", chip(value("B_gain_AFmax"), level_color(value("B_gain_AFmax"), c(0.01, 0.001, 0.0001), high_is_bad = FALSE))),
                      kv("Benign loss AF max", chip(value("B_loss_AFmax"), level_color(value("B_loss_AFmax"), c(0.01, 0.001, 0.0001), high_is_bad = FALSE)))),
                 card("Links", "external-link-alt", "#5e35b1",
                      div(class = "pr-links", lapply(seq_len(nrow(links)), function(i)
                        tags$a(href = links$url[i], target = "_blank", title = links$description[i], class = "btn btn-xs pr-link-variant", links$name[i]))))),
          column(4,
                 card(sprintf("Genes (%d)", length(genes)), "dna", "#1565c0",
                      div(style = "max-height: 110px; overflow-y: auto; word-break: break-word; font-size: 13px;", paste(genes, collapse = ", "))),
                 card("ACMG criteria (AnnotSV)", "balance-scale", "#ef6c00",
                      div(style = "font-size: 13px; word-break: break-word;", lapply(strsplit(value("AnnotSV_ranking_criteria"), ";", fixed = TRUE)[[1]], div))),
                 card("Gene constraint", "shield-alt", "#455a64",
                      kv("LOEUF bin", chip(value("LOEUF_bin"), level_color(value("LOEUF_bin"), c(5, 3, 1), high_is_bad = FALSE))),
                      kv("pLI", chip(value("GnomAD_pLI"), level_color(value("GnomAD_pLI"), c(0.5, 0.75, 0.9)))),
                      kv("Haploinsufficiency (HI)", value("HI")), kv("Triplosensitivity (TS)", value("TS")))),
          column(4,
                 card("Known pathogenic CNVs", "exclamation-triangle", "#c62828",
                      strong("Gains"), long_text("P_gain_phen"), strong("Losses"), long_text("P_loss_phen")),
                 card("OMIM / GenCC", "notes-medical", "#ad1457", long_text("OMIM_phenotype"), strong("GenCC"), long_text("GenCC_disease")))))
    })

    output$download <- downloadHandler(
      filename = function() paste0("priorr_cnvs_", Sys.Date(), ".tsv"),
      content = function(file) data.table::fwrite(filtered()[, display_columns(), drop = FALSE], file, sep = "\t", na = ""))

    filtered
  })
}
