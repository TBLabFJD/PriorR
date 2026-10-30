# SNV analysis: sidebar with the filters and profiles, variant table, variant details and download.
# By default nothing is filtered; profiles (built-in or saved by each user) set the filters and
# the displayed columns.

snv_sidebar_ui <- function(id) {
  ns <- NS(id)
  picker <- function(id, label, choices = NULL, selected = choices)
    shinyWidgets::pickerInput(ns(id), label, choices = choices, selected = selected,
                              options = list(`actions-box` = TRUE, `live-search` = TRUE, size = 10), multiple = TRUE)
  check <- function(id, label, value = FALSE)
    shinyWidgets::prettyCheckbox(ns(id), label, value = value, outline = TRUE, status = "success")
  file <- function(id, label, accept = c(".txt", ".csv", ".tsv")) fileInput(ns(id), label, accept = accept)
  panels <- function(id, label) selectizeInput(ns(id), label, choices = NULL, multiple = TRUE,
                                               options = list(placeholder = "Search a panel"))

  shinydashboard::sidebarMenu(
    id = ns("menu"),
    shinydashboard::menuItem("Load variants", icon = icon("upload"), startExpanded = TRUE,
      load_button_ui(ns("load"))),
    snv_project_sidebar_ui(ns),
    shinydashboard::menuItem("Profiles", icon = icon("sliders-h"),
      selectInput(ns("profile"), "Filtering profile", choices = profile_choices()),
      actionButton(ns("apply_profile"), "Apply", icon = icon("check"), style = "margin: 0 15px;"),
      textInput(ns("profile_name"), "Save the current filters and columns as"),
      actionButton(ns("save_profile"), "Save profile", icon = icon("save"), style = "margin: 0 15px;"),
      actionLink(ns("delete_profile"), "Delete the selected profile (only yours)", style = "margin: 6px 15px; display: block;")),
    shinydashboard::menuItem("Columns", icon = icon("columns"),
      picker("columns", NULL)),
    shinydashboard::menuItem("Filters", icon = icon("ruler-horizontal"),
      picker("type", "Type"),
      picker("region", "Region", c(genomic_regions, "NA")),
      picker("consequence", "Consequence"),
      file("bed", "Regions (BED)", c(".bed", ".txt", ".tsv")),
      check("canonical", "Canonical transcripts only"),
      check("exclude_predicted", "Exclude predicted transcripts (XM_, XR_)")),
    shinydashboard::menuItem("Frequencies", icon = icon("users"),
      sliderInput(ns("max_gnomadg_af"), "Max gnomADg AF", min = 0, max = 1, value = 1),
      sliderInput(ns("max_gnomadg_popmax"), "Max gnomADg AF popmax", min = 0, max = 1, value = 1),
      sliderInput(ns("max_gnomade_af"), "Max gnomADe AF", min = 0, max = 1, value = 1),
      sliderInput(ns("max_fjd_af"), "Max FJD_MAF AF", min = 0, max = 1, value = 1)),
    shinydashboard::menuItem("Pathogenicity", icon = icon("file-medical-alt"),
      shinyWidgets::sliderTextInput(ns("min_cadd"), "Min CADD PHRED", choices = seq(0, 75, by = 5), selected = 0),
      uiOutput(ns("score_filters")),
      radioButtons(ns("scores_mode"), "The scores above (except CADD)", c("all must pass (missing scores pass)" = "all",
                                                                          "at least one must pass" = "any")),
      check("exclude_benign", "Exclude ClinVar benign")),
    shinydashboard::menuItem("Inheritance", icon = icon("child"),
      selectInput(ns("inheritance"), "Inheritance pattern", c("(any)" = "", inheritance_modes)),
      helpText("With a PED file (Family / ROH) the whole family is used.")),
    shinydashboard::menuItem("Genes and panels", icon = icon("dna"),
      textInput(ns("genes"), "Your genes (comma separated)"),
      panels("panels", "Gene panels (PanelApp, RetiGene, Orphanet, saved)"),
      check("amber", "Include amber genes of PanelApp"),
      check("red", "Include amber and red genes of PanelApp"),
      uiOutput(ns("red_note")),
      panels("specific_panels", "Disease-specific genes (Benítez et al., 2026)"),
      helpText("Genes predicted by GLOWgenes as specific of a disease or of a disease class; they are added to the genes of the panels."),
      file("panel_file", "...or upload a panel (one gene per line)"),
      textInput(ns("panel_name"), "Save the uploaded panel in the library as"),
      actionButton(ns("save_panel"), "Save panel", icon = icon("save"), style = "margin: 0 15px;"),
      panels("excluded_panels", "Exclude the genes of these panels"),
      file("excluded_file", "...or upload genes to exclude"),
      uiOutput(ns("panelapp_info")),
      actionButton(ns("update_panelapp"), "Download the latest PanelApp (only for me)", icon = icon("sync"), style = "margin: 0 15px;")),
    shinydashboard::menuItem("GLOWgenes ranking", icon = icon("sort-amount-down"),
      radioButtons(ns("glow_source"), "Rank the genes with",
                   c("the GLOWgenes column of the file" = "column", "precomputed GLOWgenes rankings" = "panel",
                     "a new ranking (run GLOWgenes)" = "run")),
      conditionalPanel("input.glow_source == 'panel'", ns = ns,
        selectizeInput(ns("glow_panels"), "Rankings (the best rank of a gene when several)", choices = glow_rank_choices(), multiple = TRUE,
                       options = list(placeholder = "Orphanet, PanelApp or computed by the lab", maxOptions = 1000))),
      conditionalPanel("input.glow_source == 'run'", ns = ns,
        uiOutput(ns("glow_docker")),
        selectizeInput(ns("glow_run_panel"), sprintf("Genes of a panel (%d or more genes)", glow_min_genes),
                       choices = panel_choices(min_genes = glow_min_genes), options = list(placeholder = "Search a panel", maxOptions = 2000,
                                      onInitialize = I("function() { this.setValue(''); }"))),
        fileInput(ns("glow_run_file"), "...or a file of genes (one per line)", accept = c(".txt", ".csv", ".tsv")),
        textInput(ns("glow_run_name"), "Name of the ranking"),
        actionButton(ns("glow_run"), "Run GLOWgenes (Docker)", icon = icon("play"), style = "margin: 0 15px;"),
        uiOutput(ns("glow_run_status")),
        helpText("When it finishes, the new ranking is added to the precomputed rankings and selected.")),
      check("glow_limit", "Only genes up to a GLOWgenes rank"),
      conditionalPanel("input.glow_limit", ns = ns,
        sliderInput(ns("glow"), "Max GLOWgenes rank", min = 0, max = 1800, value = 100)),
      helpText("The rankings are shown with the 'GLOWgenes ranking' view (button above the table); the max rank also filters the sidebar view.")),
    shinydashboard::menuItem("Phenotypes / HPO", icon = icon("user-cog"),
      selectizeInput(ns("phenotype"), "HPO terms", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Type a phenotype or HP: identifier")),
      textInput(ns("hpo_ids"), "...or paste HPO identifiers", placeholder = "HP:0000510, HP:0000556"),
      helpText(sprintf("The genes of the terms (and of their subterms) work like a gene panel: they are added to the genes of the panels in Genes and panels. HPO %s.", hpo_version()))),
    shinydashboard::menuItem("Family / ROH", icon = icon("home"),
      radioButtons(ns("family_source"), "Family", c("none" = "none", "trio (samples of the file)" = "trio", "PED file" = "ped")),
      conditionalPanel("input.family_source == 'trio'", ns = ns,
        selectInput(ns("trio_proband"), "Proband (affected)", choices = NULL),
        radioButtons(ns("trio_sex"), "Sex of the proband", c("male" = "1", "female" = "2"), inline = TRUE),
        selectInput(ns("trio_father"), "Father", choices = NULL),
        selectInput(ns("trio_mother"), "Mother", choices = NULL),
        helpText("Choose 'De novo' (or another pattern) in Inheritance.")),
      conditionalPanel("input.family_source == 'ped'", ns = ns,
        file("ped", "Family (PED file)", c(".ped", ".txt", ".tsv"))),
      textInput(ns("proband"), "Proband (for ROH)"),
      check("roh", "Color variants in ROH")),
    shinydashboard::menuItem("VCF from other test (same patient)", icon = icon("project-diagram"),
      open_file_ui(ns("other_test"), "Open VCF of another test", "Open the VCF of another test (.vcf / .vcf.gz)"),
      check("other_test_colour", "Highlight variants NOT in the other test"),
      check("other_test_filter", "Hide variants already in the other test")),
    shinydashboard::menuItem("Analysis rounds", icon = icon("layer-group"),
      helpText("A round keeps the variants of the current gene panels / HPO / genes, frequency and CADD (the other filters are not used). With new filters, hide them to see only the new variants."),
      textInput(ns("round_name"), "Name of the round (empty: sample and filters)"),
      actionButton(ns("save_round"), "Save a round with the current panels, frequency and CADD", icon = icon("save"), style = "margin: 0 15px;"),
      check("hide_rounds", "Hide the variants of previous rounds"),
      uiOutput(ns("rounds"))),
    shinydashboard::menuItem("Classification", icon = icon("th-list"),
      check("annotate", "Show classifications and comments"),
      check("reanalysis", "Only classified variants with changes since their classification (reanalysis)")),
    shinydashboard::menuItem("Save data", icon = icon("save"),
      downloadButton(ns("download"), "Download filtered data", class = "btn-primary", style = "margin: 6px 15px;")))
}

# Patient of the loaded file, for the header of the page (output of snv_server)
snv_patient_ui <- function(id) uiOutput(NS(id)("patient"), inline = TRUE)

snv_main_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "pr-viewbar", style = "display: flex; align-items: center; gap: 15px; margin-bottom: 6px;",
        shinyWidgets::radioGroupButtons(ns("view"), NULL, status = "success", size = "sm",
          choices = c("Sidebar filters" = "filters", "PriorR priority" = "priority", "GLOWgenes ranking" = "glow"), selected = "filters"),
        conditionalPanel("input.view == 'filters'", ns = ns, style = "margin-bottom: -10px;",
          shinyWidgets::prettyCheckbox(ns("color_priority"), "Color red / orange", value = FALSE, outline = TRUE, status = "success")),
        conditionalPanel("input.view == 'glow'", ns = ns, style = "margin-bottom: -10px;",
          shinyWidgets::prettyCheckbox(ns("glow_hide_rank0"), "Hide panel genes", value = FALSE, outline = TRUE, status = "success")),
        uiOutput(ns("loaded_info"), inline = TRUE)),
    uiOutput(ns("out_info")),
    uiOutput(ns("glow_info")),
    uiOutput(ns("priority_info")),
    uiOutput(ns("reanalysis_info")),
    DT::DTOutput(ns("table")),
    br(),
    variant_info_ui(ns("info")),
    reads_ui(ns("reads")),
    load_modal_ui(ns("load")))
}

snv_server <- function(id, user) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    data <- reactiveVal(NULL)
    applied_profile <- reactiveVal(NULL)
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
    updateSelectizeInput(session, "panels", choices = panel_choices(specific = FALSE), server = TRUE)
    updateSelectizeInput(session, "specific_panels", choices = panel_choices(specific = TRUE), server = TRUE)
    updateSelectizeInput(session, "excluded_panels", choices = panel_choices(), server = TRUE)

    loader <- load_server("load", user, function() {
      p <- current_profile()
      list(filters = p$filters, columns = p$columns)
    })
    # a profile saved or deleted in the load dialog appears in the sidebar too
    observeEvent(loader$profile_saved(), {
      choices <- profile_choices(user()$id)
      updateSelectInput(session, "profile", choices = choices,
                        selected = if (isTRUE(input$profile %in% choices)) input$profile else choices[1])
    }, ignoreInit = TRUE)

    # ---- Profiles ----
    refresh_profiles <- function(selected = isolate(input$profile)) {
      choices <- profile_choices(isolate(user())$id)
      updateSelectInput(session, "profile", choices = choices, selected = if (isTRUE(selected %in% choices)) selected else choices[1])
      loader$refresh_profiles()
    }
    observeEvent(user(), refresh_profiles())

    all_or <- function(x, choices) x %||% choices
    # Score filters of the loaded file (score_filters.R); a profile sets their values
    score_preset <- reactiveVal(list())
    score_input <- function(col) paste0("score_", make.names(col))
    output$score_filters <- renderUI({
      req(data())
      defs <- score_filters_for(names(data()))
      preset <- score_preset()
      lapply(seq_len(nrow(defs)), function(i) {
        d <- defs[i, ]
        sliderInput(ns(score_input(d$col)), d$label, min = d$min, max = d$max, step = d$step,
                    value = preset[[d$col]] %||% score_off_value(d))
      })
    })
    active_scores <- reactive({
      req(data())
      defs <- score_filters_for(names(data()))
      out <- lapply(seq_len(nrow(defs)), function(i) {
        d <- defs[i, ]
        value <- input[[score_input(d$col)]]
        if (is.null(value) || value == score_off_value(d)) NULL else list(col = d$col, op = d$op, value = value)
      })
      Filter(Negate(is.null), out)
    })
    # Sets the sidebar filters and the columns of a profile (columns, types... need the data)
    apply_profile <- function(p, df = isolate(data())) {
      f <- p$filters
      shinyWidgets::updatePrettyCheckbox(session, "color_priority", value = isTRUE(f$color_priority))
      if (isTRUE(f$prioritize)) shinyWidgets::updateRadioGroupButtons(session, "view", selected = "priority")
      for (i in c("canonical", "exclude_predicted", "exclude_benign", "amber", "red"))
        shinyWidgets::updatePrettyCheckbox(session, i, value = isTRUE(f[[i]]))
      for (i in c("max_gnomadg_af", "max_gnomadg_popmax", "max_gnomade_af", "max_fjd_af"))
        updateSliderInput(session, i, value = f[[i]])
      # profiles keep glow = -1 for "no max rank"
      glow_limit <- !is.null(f$glow) && f$glow >= 0
      shinyWidgets::updatePrettyCheckbox(session, "glow_limit", value = glow_limit)
      if (glow_limit) updateSliderInput(session, "glow", value = f$glow)
      shinyWidgets::updateSliderTextInput(session, "min_cadd", selected = f$min_cadd)
      score_preset(f$scores %||% list())
      updateRadioButtons(session, "scores_mode", selected = f$scores_mode %||% "all")
      updateSelectInput(session, "inheritance", selected = f$inheritance %||% "")
      # profiles keep the specific genes panels with the other panels
      chosen <- unlist(f$panels)
      specific <- startsWith(chosen %||% character(0), "glowspecific:")
      updateSelectizeInput(session, "panels", choices = panel_choices(specific = FALSE), selected = chosen[!specific], server = TRUE)
      updateSelectizeInput(session, "specific_panels", choices = panel_choices(specific = TRUE), selected = chosen[specific], server = TRUE)
      updateSelectizeInput(session, "excluded_panels", choices = panel_choices(), selected = unlist(f$excluded_panels), server = TRUE)
      updateSelectizeInput(session, "phenotype", choices = hpo_choices(), selected = unlist(f$phenotypes) %||% character(0), server = TRUE)
      shinyWidgets::updatePickerInput(session, "region", selected = all_or(unlist(f$region), c(genomic_regions, "NA")))
      if (!is.null(df)) {
        consequences <- c(sort(unique(df$Consequence)), "NA")
        types <- sort(unique(df$VARIANT_CLASS))
        shinyWidgets::updatePickerInput(session, "consequence", choices = consequences, selected = all_or(unlist(f$consequence), consequences))
        shinyWidgets::updatePickerInput(session, "type", choices = types, selected = all_or(unlist(f$type), types))
        sample_cols <- grep("_(GT|VAF|AD|DP)$", names(df), value = TRUE)
        default_columns <- intersect(c(default_snv_columns[1:8], sample_cols, default_snv_columns[-(1:8)]), names(df))
        columns <- intersect(unlist(p$columns), names(df))
        if (!length(columns)) columns <- default_columns
        shinyWidgets::updatePickerInput(session, "columns", choices = ordered_column_choices(names(df), columns),
                                        selected = columns)
      }
      applied_profile(p$name)
    }

    observeEvent(input$apply_profile, {
      p <- all_profiles(user()$id)[[input$profile]]
      req(p)
      apply_profile(p)
      showNotification(paste("Profile", p$name, "applied"), type = "message")
    })

    # Current filters and columns as a profile ("all selected" is saved as no filter, so the
    # profile works with any file)
    current_profile <- function() {
      df <- data()
      if_not_all <- function(x, choices) if (length(choices) && setequal(x, choices)) NULL else as.list(x)
      list(
        name = input$profile_name,
        description = paste("Saved by", user()$username, "on", Sys.Date()),
        load = loader$settings(),
        filters = list(
          type = if_not_all(input$type, if (is.null(df)) NULL else unique(df$VARIANT_CLASS)),
          region = if_not_all(input$region, c(genomic_regions, "NA")),
          consequence = if_not_all(input$consequence, if (is.null(df)) NULL else c(unique(df$Consequence), "NA")),
          canonical = input$canonical, exclude_predicted = input$exclude_predicted,
          max_gnomadg_af = input$max_gnomadg_af, max_gnomadg_popmax = input$max_gnomadg_popmax,
          max_gnomade_af = input$max_gnomade_af, max_fjd_af = input$max_fjd_af,
          min_cadd = as.numeric(input$min_cadd), exclude_benign = input$exclude_benign,
          scores = if (!is.null(df)) stats::setNames(lapply(active_scores(), `[[`, "value"), vapply(active_scores(), `[[`, "", "col")),
          scores_mode = input$scores_mode,
          inheritance = input$inheritance, panels = as.list(c(input$panels, input$specific_panels)), amber = input$amber, red = isTRUE(input$red),
          excluded_panels = as.list(input$excluded_panels), phenotypes = as.list(input$phenotype), glow = if (isTRUE(input$glow_limit)) input$glow else -1,
          color_priority = isTRUE(input$color_priority) || identical(input$view, "priority"), prioritize = identical(input$view, "priority")),
        columns = as.list(input$columns))
    }

    observeEvent(input$save_profile, {
      result <- tryCatch(save_profile(user()$id, current_profile()), error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      refresh_profiles(paste0("user:", result))
      updateTextInput(session, "profile_name", value = "")
      showNotification(paste("Profile", result, "saved"), type = "message")
    })

    observeEvent(input$delete_profile, {
      if (!startsWith(input$profile, "user:")) return(showNotification("Built-in profiles can't be deleted", type = "warning"))
      delete_profile(user()$id, sub("^user:", "", input$profile))
      refresh_profiles()
    })

    # ---- Data ----
    # Selections of a restored session, applied when its data is loaded again
    restored <- reactiveValues(input = NULL)
    onRestore(function(state) {
      restored$input <- state$input[c("columns", "consequence", "type", "region", "min_cadd")]
    })

    observeEvent(loader$loaded(), {
      l <- loader$loaded()
      df <- add_annotations(l$df)
      data(df)
      r <- restored$input
      restored$input <- NULL
      if (isTRUE(l$restoring) && !is.null(r)) {
        shinyWidgets::updatePickerInput(session, "columns", choices = ordered_column_choices(names(df), r$columns), selected = r$columns)
        shinyWidgets::updatePickerInput(session, "consequence", choices = c(sort(unique(df$Consequence)), "NA"), selected = r$consequence)
        shinyWidgets::updatePickerInput(session, "type", choices = sort(unique(df$VARIANT_CLASS)), selected = r$type)
        shinyWidgets::updatePickerInput(session, "region", selected = r$region)
        shinyWidgets::updateSliderTextInput(session, "min_cadd", selected = r$min_cadd)
      } else {
        apply_profile(l$profile, df)
      }
    })

    # Patient of the loaded file, in the header of the page: the proband (the first sample, or the
    # one chosen for the trio), the other samples and the file
    output$patient <- renderUI({
      req(data())
      samples <- sample_names(names(data()))
      proband <- if (isTRUE(input$trio_proband %in% samples)) input$trio_proband else samples[1]
      file <- loader$loaded()$file
      patient_header(proband, setdiff(samples, proband), if (length(file)) basename(file))
    })

    output$loaded_info <- renderUI({
      req(data(), filtered())
      hidden <- length(hidden_keys())
      n <- function(x) format(x, big.mark = ",")
      out <- showing_out()
      # the two sets of variants: those that pass the filters and those that don't (links to show each one)
      set_link <- function(id, label, active, color)
        actionLink(ns(id), label, class = "pr-set-link", style = sprintf("color: %s;%s", color, if (active) " font-weight: 700; text-decoration: underline;" else ""))
      # rows of the file, rows that passed the load filters (all loaded), and then the sidebar filters
      total <- loader$loaded()$total_rows
      loaded_text <- if (length(total) && total > nrow(data()))
        sprintf("%s rows in the file · %s passed the load filters · ", n(total), n(nrow(data())))
      else sprintf("%s rows loaded · ", n(nrow(data())))
      p(class = "text-muted", style = "margin: 0;", title = "pass / filtered out: the sidebar filters (or the view) on the loaded rows",
        loaded_text,
        set_link("show_pass", paste(n(nrow(filtered())), "pass"), !out, "#2e7d32"), " · ",
        set_link("show_out", paste(n(length(out_rows())), "filtered out"), out, "#c62828"),
        if (!is.null(project())) sprintf(" · project %s: %s other samples", project()$name, n(project()$total)),
        sprintf(" · profile: %s%s", applied_profile() %||% "None",
                if (hidden) sprintf(" · hiding %s variants of previous rounds", n(hidden)) else ""))
    })

    # ---- Variants that don't pass the filters: shown instead of the others until going back ----
    showing_out <- reactiveVal(FALSE)
    out_rows <- reactive(excluded_rows(data(), filtered(), by_variant = identical(input$view, "glow")))
    observeEvent(input$show_out, showing_out(TRUE))
    observeEvent(input$show_pass, showing_out(FALSE))
    observeEvent(input$back_to_pass, showing_out(FALSE))
    observeEvent(data(), showing_out(FALSE))
    # the table of the page: the variants that pass, or those that don't
    table_rows <- reactive({
      if (!showing_out()) return(filtered())
      excluded_snv(data(), out_rows(), params(), kept = filtered())
    })
    output$out_info <- renderUI({
      req(showing_out(), data())
      div(style = "padding: 8px 14px; margin-bottom: 8px; background: #fdecea; color: #b71c1c; border-left: 4px solid #e53935; border-radius: 4px;",
          icon("filter"),
          sprintf(" Showing the %s variants that don't pass the filters. ", format(length(out_rows()), big.mark = ",")),
          actionLink(ns("back_to_pass"), sprintf("Back to the %s that pass", format(nrow(filtered()), big.mark = ",")),
                     style = "font-weight: 700; color: #1b5e20;"))
    })

    # ---- Analysis rounds (snv_rounds.R) ----
    hidden_keys <- snv_rounds_server(input, output, session, user, data, function() params(), loader$loaded)

    # ---- Project: counts in the other samples (snv_project.R) ----
    project <- snv_project_server(input, output, session, user, data, loader$loaded, loader$open_path)

    # ---- PanelApp of the user (snv_panelapp.R) ----
    snv_panelapp_server(input, output, session, user)

    # ---- Side files ----
    uploaded_panel <- reactive(if (is.null(input$panel_file)) NULL else read_gene_list(input$panel_file$datapath))
    observeEvent(input$save_panel, {
      result <- tryCatch(save_custom_panel(input$panel_name, uploaded_panel()), error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      updateSelectizeInput(session, "panels", choices = panel_choices(specific = FALSE), selected = isolate(input$panels), server = TRUE)
      updateSelectizeInput(session, "excluded_panels", choices = panel_choices(), selected = isolate(input$excluded_panels), server = TRUE)
      updateSelectizeInput(session, "glow_run_panel", choices = panel_choices(min_genes = glow_min_genes), selected = isolate(input$glow_run_panel))
      updateTextInput(session, "panel_name", value = "")
      showNotification(paste("Panel", input$panel_name, "saved in the library"), type = "message")
    })
    bed <- reactive(if (is.null(input$bed)) NULL else read_bed(input$bed$datapath))
    # Family: a trio chosen among the samples of the file, or a PED file
    ped <- reactive({
      source <- input$family_source %||% (if (is.null(input$ped)) "none" else "ped")
      if (source == "trio") {
        if (!nzchar(input$trio_proband %||% "")) return(NULL)
        trio_ped(input$trio_proband, input$trio_father, input$trio_mother, input$trio_sex %||% "1")
      } else if (source == "ped" && !is.null(input$ped)) read_ped(input$ped$datapath)
      else NULL
    })
    observeEvent(input$ped, updateRadioButtons(session, "family_source", selected = "ped"))
    # De novo needs both parents
    observeEvent(list(input$inheritance, ped()), {
      req(identical(input$inheritance, "De novo"), data())
      fam <- if (is.null(ped())) NULL else ped_samples(ped(), names(data()))
      if (is.null(fam) || is.null(fam$father) || is.null(fam$mother) || !length(fam$affected))
        showNotification("De novo needs the proband and both parents: choose the trio in Family / ROH (or a PED file).",
                         type = "warning", id = "denovo_trio")
    }, ignoreInit = TRUE)
    other_test_file <- open_file_server("other_test", input, output)
    other_test_ids <- reactive(if (is.null(other_test_file())) NULL else read_vcf_ids(other_test_file()))

    # Classifications and comments visible to the user
    info <- NULL
    annotations <- reactive({
      req(isTRUE(input$annotate), user())
      info$changed()
      comment_annotations(user()$id)
    })

    # Reanalysis: classified variants of this file whose annotation changed since they were classified
    reanalysis_table <- reactive({
      req(data(), user())
      info$changed()
      reanalysis(data(), visible_comments(user()$id))
    })
    output$reanalysis_info <- renderUI({
      r <- reanalysis_table()
      req(nrow(r) > 0)
      div(class = "alert alert-info", style = "padding: 8px 12px;",
          icon("redo"), strong(sprintf(" Reanalysis: %d classified variant%s changed since %s classification. ", nrow(r),
                                       if (nrow(r) > 1) "s have" else " has", if (nrow(r) > 1) "their" else "its")),
          if (!isTRUE(input$reanalysis)) "Tick 'Only classified variants with changes' in Classification to see them.",
          tags$details(tags$summary("Changes"),
                       tags$ul(style = "margin: 4px 0 0 0;", lapply(seq_len(nrow(r)), function(i)
                         tags$li(strong(r$gene[i]), sprintf("%s (%s, %s by %s): %s", r$variant_key[i], r$classification[i],
                                                            r$date[i], r$author[i], r$changes[i]))))))
    })

    # ---- Filtering ----
    params <- reactive({
      roh_samples <- if (nzchar(input$proband)) input$proband
                     else if (!is.null(ped())) ped()$V2[ped()$V6 == "2"]
      # the genes of the HPO terms work like one more panel
      panel <- c(panel_genes(c(input$panels, input$specific_panels), amber = isTRUE(input$amber), red = isTRUE(input$red)), uploaded_panel(),
                 hpo_genes(input$phenotype))
      excluded <- c(panel_genes(input$excluded_panels, amber = TRUE),
                    if (!is.null(input$excluded_file)) read_gene_list(input$excluded_file$datapath))
      snv_params(
        consequence = input$consequence, region = input$region, type = input$type,
        canonical = input$canonical, exclude_predicted = input$exclude_predicted,
        genes = split_genes(input$genes),
        max_af = list(gnomADg_AF = input$max_gnomadg_af, gnomADg_AF_popmax = input$max_gnomadg_popmax,
                      gnomADe_AF = input$max_gnomade_af, FJD_MAF_AF = input$max_fjd_af),
        min_cadd = as.numeric(input$min_cadd),
        scores = if (is.null(data())) NULL else active_scores(), scores_mode = input$scores_mode %||% "all",
        inheritance = input$inheritance, ped = ped(),
        glow_max = if (isTRUE(input$glow_limit)) input$glow, glow_sources = glow_panels(),
        panel_genes = if (length(panel)) unique(panel), excluded_genes = if (length(excluded)) unique(excluded),
        bed = bed(),
        other_test_ids = other_test_ids(), hide_other_test = input$other_test_filter,
        only_keys = if (isTRUE(input$reanalysis)) reanalysis_table()$variant_key,
        hidden_keys = hidden_keys(),
        exclude_clinvar_benign = input$exclude_benign, prioritize = identical(input$view, "priority"),
        annotations = if (isTRUE(input$annotate)) annotations(),
        columns = input$columns, roh_samples = roh_samples,
        project = project(), max_project = if (!is.null(project()) && isTRUE(input$max_project >= 0)) input$max_project)
    })

    # GLOWgenes ranking view: its own pathogenicity filters (glowgenes.json) instead of the sidebar ones
    glow_panels <- reactive({
      if (!identical(input$glow_source, "column") || !has_glow_column(data())) input$glow_panels else NULL
    })

    # ---- Run GLOWgenes (snv_glow_run.R) ----
    snv_glow_run_server(input, output, session, user)
    glow_ready <- reactive(has_glow_column(data()) && identical(input$glow_source, "column") || length(input$glow_panels) > 0)

    filtered <- reactive({
      req(data(), user())
      if (identical(input$view, "glow")) {
        req(glow_ready())
        p <- params()
        return(filter_snv(data(), snv_params(columns = p$columns, annotations = p$annotations, other_test_ids = p$other_test_ids,
                                             project = p$project, max_project = p$max_project,
                                             hide_other_test = p$hide_other_test, roh_samples = p$roh_samples, glow_max = p$glow_max, hidden_keys = p$hidden_keys,
                                             glow = list(sources = glow_panels(), hide_rank0 = isTRUE(input$glow_hide_rank0)))))
      }
      filter_snv(data(), params())
    })

    # ---- Boxes of the PriorR priority and GLOWgenes views (snv_view_info.R) ----
    snv_view_info_server(input, output, session, data, filtered, glow_ready, glow_panels)

    # ---- Table (processed in the server: the browser only gets the visible page) ----
    # columns of the table shown, to keep its filter boxes, search and order when it is drawn again
    # (other columns, filters or view); a new file starts clean
    table_columns <- reactiveVal(NULL)
    observeEvent(data(), table_columns(NULL), priority = 10)
    output$table <- DT::renderDT({
      DF <- table_rows()
      shown <- display_table(DF)
      shown$Classify <- classify_buttons(nrow(shown), ns("classify"))
      colors <- row_colors(DF, roh = input$roh, other_test = input$other_test_colour,
                           priority = isTRUE(input$color_priority) || identical(input$view, "priority") && !showing_out())
      state <- isolate(table_view_state(table_columns(), input$table_search_columns, input$table_search, input$table_view))
      table_columns(c(names(shown), ".row_color"))
      priorr_datatable(shown, ns("full_cell"), colors, reorder_input = ns("column_order"), state = state)
    }, server = TRUE, funcFilter = priorr_dt_filter)
    # Dragging a header moves the column in the picker, so the table, the download and the
    # profiles keep the new order
    observeEvent(input$column_order, {
      df <- data()
      req(df)
      columns <- reordered_columns(input$columns %||% names(df), input$column_order)
      shinyWidgets::updatePickerInput(session, "columns", choices = ordered_column_choices(names(df), columns), selected = columns)
    })
    observeEvent(input$full_cell, show_full_cell(input$full_cell))

    selected <- reactive({
      i <- input$table_rows_selected
      if (length(i) != 1 || i > nrow(table_rows())) NULL else table_rows()[i, , drop = FALSE]
    })
    info <- variant_info_server("info", selected, data, user, ped, project, table_rows, reactive(list(phenotype = input$phenotype)))
    reads_server("reads", selected, data)

    # "Classify" on a row: select it and open its comments
    observeEvent(input$classify, {
      i <- input$classify
      DT::selectRows(DT::dataTableProxy("table"), i)
      v <- table_rows()[i, ]
      info$open(list(key = variant_key(v$CHROM, v$POS, v$REF, v$ALT), gene = as.character(v$SYMBOL), label = variant_label(v)))
    })

    output$download <- downloadHandler(
      filename = function() paste0("priorr_variants_", Sys.Date(), ".tsv"),
      content = function(file) data.table::fwrite(display_table(table_rows()), file, sep = "\t", na = ""))

    setBookmarkExclude(c("apply_profile", "save_profile", "delete_profile", "profile_name", "save_panel", "panel_name"))
    observeEvent(data(), {
      # without GLOWgenes column the ranking needs an Orphanet panel
      if (!has_glow_column(data())) updateRadioButtons(session, "glow_source", selected = "panel")
      # samples for the trio
      samples <- sample_names(names(data()))
      updateSelectInput(session, "trio_proband", choices = samples, selected = samples[1])
      updateSelectInput(session, "trio_father", choices = c("(none)" = "", samples), selected = if (length(samples) >= 2) samples[2] else "")
      updateSelectInput(session, "trio_mother", choices = c("(none)" = "", samples), selected = if (length(samples) >= 3) samples[3] else "")
    })

    filtered
  })
}
