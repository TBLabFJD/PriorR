# Loading of a variant file (gene panel, exome or genome: always the same way).
# The file is read with the WGS prefilter (wgs_prefilter.R) in a background R process: the
# load filters of the chosen profile keep only the candidates (none by default), and the big
# files are cached so that loading them again takes seconds.
# Returns list(loaded = reactive list(df, profile), settings = reactive load filters,
#              refresh_profiles = function()).

load_button_ui <- function(id) {
  shinyWidgets::actionBttn(NS(id, "open"), "Open variant file", icon = icon("folder-open"),
                           style = "minimal", color = "success")
}

load_modal_ui <- function(id) {
  ns <- NS(id)
  check <- function(id, label) shinyWidgets::prettyCheckbox(ns(id), label, outline = TRUE, status = "success")
  catalog_choices <- panel_choices()
  static_modal(
    ns("dialog"), "Open a variant file",
    fluidRow(
      column(6,
             h4("1. File"),
             open_file_ui(ns("file"), "Open file", "Open a variant file (.tsv / .tsv.gz): panel, exome or genome", sidebar = FALSE),
             h4("2. Filtering profile"),
             radioButtons(ns("profile"), NULL, choices = profile_choices()),
             uiOutput(ns("profile_description")),
             div(style = "background: #f4f6f8; padding: 8px 10px; border-radius: 4px; margin-bottom: 10px;",
                 fluidRow(
                   column(8, textInput(ns("profile_name"), "Save the current filters as a profile", placeholder = "Name")),
                   column(4, br(), actionButton(ns("save_profile"), "Save", icon = icon("save")))),
                 fluidRow(
                   column(8, fileInput(ns("profile_file"), "Upload a profile (.json)", accept = ".json")),
                   column(4, br(), downloadButton(ns("download_profile"), "Download"))),
                 actionLink(ns("delete_profile"), "Delete the selected profile (only yours)")),
             h4("3. Genes and panels"),
             textInput(ns("genes"), "Only these genes (comma separated)", placeholder = "All genes"),
             selectizeInput(ns("panels"), "...and the genes of these panels", choices = catalog_choices, multiple = TRUE,
                            options = list(placeholder = "Click to see the panels (PanelApp, Orphanet, saved)", maxOptions = 2000)),
             check("amber", "Include PanelApp amber genes"),
             check("red", "Include PanelApp amber and red genes"),
             fileInput(ns("panel_file"), "...and the genes of a panel file (one gene per line)", accept = c(".txt", ".csv", ".tsv")),
             selectizeInput(ns("excluded_panels"), "Exclude the genes of these panels", choices = catalog_choices, multiple = TRUE,
                            options = list(placeholder = "None", maxOptions = 2000))),
      column(6,
             h4("4. Load filters"),
             helpText("Only the variants that pass these filters are loaded. Empty = no filter."),
             shinyWidgets::pickerInput(ns("regions"), "Genomic regions", choices = genomic_regions, multiple = TRUE,
                                       options = list(`actions-box` = TRUE, `none-selected-text` = "All")),
             check("canonical", "Canonical transcripts only"),
             fluidRow(
               column(6, numericInput(ns("af_gnomadg_popmax"), "Max gnomADg AF popmax", NA, min = 0, max = 1, step = 0.001)),
               column(6, numericInput(ns("af_gnomadg"), "Max gnomADg AF", NA, min = 0, max = 1, step = 0.001))),
             fluidRow(
               column(6, numericInput(ns("af_gnomade"), "Max gnomADe AF", NA, min = 0, max = 1, step = 0.001)),
               column(6, numericInput(ns("af_fjd"), "Max FJD_MAF AF", NA, min = 0, max = 1, step = 0.001))),
             fluidRow(
               column(6, numericInput(ns("min_cadd"), "Min CADD PHRED", NA, min = 0, max = 75, step = 1)),
               column(6, br(), check("excl_benign", "Exclude ClinVar benign"))),
             strong("Always keep, even outside the filters:"),
             fluidRow(
               column(6, check("rescue_clinvar", "ClinVar pathogenic")),
               column(6, numericInput(ns("rescue_spliceai"), "SpliceAI DS max >=", NA, min = 0, max = 1, step = 0.1))))),
    hr(),
    actionButton(ns("run"), "Load", icon = icon("play"), class = "btn-success btn-lg"),
    br(), br(),
    uiOutput(ns("status")))
}

# current_view: function returning list(filters, columns) of the sidebar, saved with the profiles
load_server <- function(id, user, current_view = function() list()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    state <- reactiveValues(running = FALSE, workdir = NULL, status = NULL, profile = NULL)
    # the PanelApp copy of the user (panels.R)
    observeEvent(user(), for (i in c("panels", "excluded_panels"))
      updateSelectizeInput(session, i, choices = panel_choices(), selected = isolate(input[[i]])))
    loaded <- reactiveVal(NULL)
    restoring <- reactiveVal(FALSE)

    observeEvent(input$open, toggle_modal(ns("dialog")))
    path <- open_file_server("file", input, output)

    # ---- Profiles ----
    refresh_profiles <- function(selected = isolate(input$profile)) {
      choices <- profile_choices(isolate(user())$id)
      updateRadioButtons(session, "profile", choices = choices,
                         selected = if (isTRUE(selected %in% choices)) selected else choices[1])
    }
    observeEvent(user(), refresh_profiles())

    current_profile <- reactive({
      req(input$profile)
      all_profiles(user()$id)[[input$profile]]
    })
    output$profile_description <- renderUI(helpText(current_profile()$description))

    set_load_inputs <- function(l) {
      af <- function(col) json_num(l$max_af[[col]])
      shinyWidgets::updatePickerInput(session, "regions", selected = unlist(l$regions) %||% character(0))
      shinyWidgets::updatePrettyCheckbox(session, "canonical", value = isTRUE(l$canonical_only))
      shinyWidgets::updatePrettyCheckbox(session, "excl_benign", value = isTRUE(l$exclude_clinvar_benign))
      shinyWidgets::updatePrettyCheckbox(session, "rescue_clinvar", value = isTRUE(l$rescue_clinvar_pathogenic))
      values <- list(af_gnomadg_popmax = af("gnomADg_AF_popmax"), af_gnomadg = af("gnomADg_AF"), af_gnomade = af("gnomADe_AF"),
                     af_fjd = af("FJD_MAF_AF"), min_cadd = json_num(l$min_cadd), rescue_spliceai = json_num(l$rescue_spliceai_min))
      for (i in names(values)) updateNumericInput(session, i, value = values[[i]])
      updateTextInput(session, "genes", value = paste(unlist(l$genes), collapse = ", "))
      shinyWidgets::updatePrettyCheckbox(session, "amber", value = isTRUE(l$amber))
      shinyWidgets::updatePrettyCheckbox(session, "red", value = isTRUE(l$red))
      for (i in c("panels", "excluded_panels")) updateSelectizeInput(session, i, selected = unlist(l[[i]]) %||% character(0))
    }
    # Choosing a profile fills the load filters (not when restoring a saved session: its own values stay)
    observeEvent(input$profile, if (!restoring()) set_load_inputs(current_profile()$load), ignoreInit = TRUE)

    # Load filters of the dialog, in the format of a profile
    settings <- reactive({
      na_null <- function(x) if (is.null(x) || is.na(x)) NULL else x
      list(regions = as.list(input$regions %||% character(0)), canonical_only = isTRUE(input$canonical),
           max_af = Filter(Negate(is.null), list(gnomADg_AF_popmax = na_null(input$af_gnomadg_popmax),
                                                gnomADg_AF = na_null(input$af_gnomadg), gnomADe_AF = na_null(input$af_gnomade),
                                                FJD_MAF_AF = na_null(input$af_fjd))),
           min_cadd = na_null(input$min_cadd), exclude_clinvar_benign = isTRUE(input$excl_benign),
           rescue_clinvar_pathogenic = isTRUE(input$rescue_clinvar), rescue_spliceai_min = na_null(input$rescue_spliceai),
           genes = as.list(split_genes(input$genes)), panels = as.list(input$panels %||% character(0)),
           amber = isTRUE(input$amber), red = isTRUE(input$red), excluded_panels = as.list(input$excluded_panels %||% character(0)))
    })

    # The selected profile with the load filters of the dialog and the current sidebar filters
    profile_to_save <- function(name) {
      view <- current_view()
      list(name = name, description = paste("Saved by", user()$username, "on", Sys.Date()),
           load = settings(), filters = view$filters %||% current_profile()$filters,
           columns = view$columns %||% current_profile()$columns)
    }
    profile_saved <- reactiveVal(0)
    after_save <- function(name) {
      refresh_profiles(paste0("user:", name))
      profile_saved(profile_saved() + 1)
      showNotification(paste("Profile", name, "saved"), type = "message")
    }

    observeEvent(input$save_profile, {
      result <- tryCatch(save_profile(user()$id, profile_to_save(input$profile_name)), error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      updateTextInput(session, "profile_name", value = "")
      after_save(result)
    })

    # A profile file (e.g. downloaded by a colleague) is saved as a profile of the user
    observeEvent(input$profile_file, {
      result <- tryCatch({
        p <- jsonlite::fromJSON(input$profile_file$datapath, simplifyVector = FALSE)
        if (!is.list(p) || is.null(p$load) && is.null(p$filters)) stop("This is not a PriorR profile")
        p$name <- p$name %||% tools::file_path_sans_ext(input$profile_file$name)
        if (p$name %in% names(builtin_profiles())) p$name <- paste(p$name, "(uploaded)")
        save_profile(user()$id, complete_profile(p))
      }, error = function(e) e)
      if (inherits(result, "error")) return(showNotification(conditionMessage(result), type = "error"))
      after_save(result)
    })

    output$download_profile <- downloadHandler(
      filename = function() paste0("priorr_profile_", gsub("[^A-Za-z0-9_-]+", "_", current_profile()$name), ".json"),
      content = function(file) jsonlite::write_json(current_profile(), file, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA))

    observeEvent(input$delete_profile, {
      if (!startsWith(input$profile, "user:")) return(showNotification("Built-in profiles can't be deleted", type = "warning"))
      delete_profile(user()$id, sub("^user:", "", input$profile))
      refresh_profiles()
      profile_saved(profile_saved() + 1)
    })

    # ---- Loading ----
    start <- function() {
      if (isolate(state$running)) return()
      file <- isolate(path())
      state$file <- file
      state$workdir <- file.path(tempdir(), paste0("load_", session$token))
      dir.create(state$workdir, showWarnings = FALSE, recursive = TRUE)
      cfg <- isolate(settings())
      cfg$max_rows <- max_rows_setting   # safety limit of the computer (config.R), not a filter
      # genes typed, of the panels and of the uploaded panel; and genes to exclude
      uploaded <- if (is.null(isolate(input$panel_file))) NULL else read_gene_list(isolate(input$panel_file)$datapath)
      cfg$genes <- as.list(unique(c(unlist(cfg$genes), panel_genes(unlist(cfg$panels), amber = cfg$amber, red = isTRUE(cfg$red)), uploaded)))
      cfg$exclude_genes <- as.list(panel_genes(unlist(cfg$excluded_panels), amber = TRUE) %||% character(0))
      profile_file <- file.path(state$workdir, "profile.json")
      jsonlite::write_json(cfg, profile_file, auto_unbox = TRUE, null = "null", digits = NA)
      unlink(file.path(state$workdir, c("progress.json", "variants.tsv")))
      log <- file.path(state$workdir, "prefilter.log")
      system2(file.path(R.home("bin"), "Rscript"),
              c(shQuote(normalizePath(wgs_script)), shQuote(file), shQuote(profile_file),
                shQuote(file.path(state$workdir, "variants.tsv")), shQuote(wgs_cache_dir),
                shQuote(file.path(state$workdir, "progress.json")), wgs_memory_limit, wgs_threads),
              stdout = log, stderr = log, wait = FALSE)
      state$profile <- isolate(current_profile())
      state$status <- list(stage = "start", message = "Starting...")
      state$running <- TRUE
    }

    observeEvent(input$run, {
      req(user())
      if (is.null(path())) return(shinyalert::shinyalert("No file", "Open a variant file first.", type = "warning"))
      start()
    })

    observe({
      req(state$running)
      invalidateLater(700)
      progress_file <- file.path(state$workdir, "progress.json")
      if (!file.exists(progress_file)) return()
      status <- tryCatch(jsonlite::fromJSON(progress_file), error = function(e) NULL)
      if (is.null(status)) return()
      state$status <- status
      if (status$stage == "done") {
        state$running <- FALSE
        df <- withProgress(message = "Loading the variants...", value = 0.5,
                           read_variants(file.path(state$workdir, "variants.tsv")))
        loaded(list(df = df, profile = state$profile, restoring = isolate(restoring()), file = state$file,
                    total_rows = status$total_rows))
        restoring(FALSE)
        toggle_modal(ns("dialog"), show = FALSE)
      } else if (status$stage == "error") {
        state$running <- FALSE
      }
    })

    output$status <- renderUI({
      s <- state$status
      req(s)
      if (s$stage == "done") {
        div(style = "color: #2e7d32;", strong("Loaded. "),
            sprintf("%s of %s rows passed the load filters (%s variants) in %s s.",
                    format(s$kept_rows, big.mark = ","), format(s$total_rows, big.mark = ","),
                    format(s$variants, big.mark = ","), s$seconds), " All of them are loaded.")
      } else if (s$stage == "error") {
        div(class = "text-danger", strong("Error: "), s$message)
      } else {
        pct <- if (is.null(s$pct)) NULL else round(s$pct)
        div(icon("spinner", class = "fa-spin"), s$message,
            if (!is.null(pct)) div(class = "progress", div(class = "progress-bar progress-bar-success",
                                                           style = sprintf("width: %s%%;", pct), paste0(pct, "%"))))
      }
    })

    # Saved sessions: the file and load filters are restored and the file is loaded again
    onRestore(function(state) restoring(isTRUE(state$input$run > 0)))
    onRestored(function(state) {
      isolate(if (isTRUE(input$run > 0) && !is.null(path())) start())
    })

    setBookmarkExclude(c("save_profile", "profile_name", "profile_file", "delete_profile"))

    list(loaded = loaded, settings = settings, refresh_profiles = refresh_profiles, profile_saved = profile_saved)
  })
}
