# Projects in the SNV tab (projects.R): choose the project of the loaded sample, manage the
# projects (files, background index) and count each variant in the other samples of the project.
# Runs inside snv_server; open_path(file) opens the load dialog with a file. Returns a reactive with
# list(id, name, counts, total, exclude) or NULL.

snv_project_sidebar_ui <- function(ns) {
  shinydashboard::menuItem("Project (other samples)", icon = icon("users-viewfinder"),
    selectInput(ns("project"), "Project", choices = c("(none)" = "", project_choices())),
    conditionalPanel("input.project != ''", ns = ns,
      selectInput(ns("project_sample"), "Open a sample of the project", choices = NULL),
      actionButton(ns("open_project_sample"), "Open", icon = icon("folder-open"), style = "margin: 0 15px 6px 15px;")),
    uiOutput(ns("project_info")),
    numericInput(ns("max_project"), "Hide variants in more than this number of other samples of the project", value = NA, min = 0, step = 1),
    actionLink(ns("manage_projects"), "Manage projects (files, index)", icon = icon("folder-tree"),
               style = "margin: 6px 15px; display: block;"))
}

snv_project_server <- function(input, output, session, user, data, loaded, open_path = function(file) NULL) {
  ns <- session$ns
  changed <- reactiveVal(0)   # projects or files added / removed
  refresh_choices <- function(selected = isolate(input$project)) {
    changed(isolate(changed()) + 1)
    updateSelectInput(session, "project", choices = c("(none)" = "", project_choices()), selected = selected)
  }

  # state of a project (files, index), read again every few seconds while it is being indexed
  state_of <- function(id) {
    s <- list(status = project_status(id), indexing = project_indexing(id))
    # while indexing, or files waiting (the indexing is starting)
    if (identical(s$indexing$stage, "running") || any(!s$status$indexed & s$status$exists)) invalidateLater(3000)
    s
  }
  project_state <- reactive({ changed(); req(nzchar(input$project %||% "")); state_of(input$project) })
  # the counts are only read again when the indexed files change
  indexed_files <- reactiveVal("")
  observe({
    s <- project_state()
    sig <- paste(s$status$path[s$status$indexed], collapse = "|")
    if (!identical(sig, isolate(indexed_files()))) indexed_files(sig)
  })

  # The samples of the project, to open one of them (the load dialog, with the file chosen)
  observe({
    s <- project_state()
    st <- s$status[s$status$exists, , drop = FALSE]
    labels <- ifelse(nzchar(st$samples), paste0(st$samples, " (", basename(st$path), ")"), basename(st$path))
    current <- if (!is.null(loaded()$file)) normal_path(loaded()$file)
    updateSelectInput(session, "project_sample", choices = stats::setNames(st$path, labels),
                      selected = if (isTRUE(isolate(input$project_sample) %in% st$path)) isolate(input$project_sample)
                                 else if (isTRUE(current %in% st$path)) current)
  })
  observeEvent(input$open_project_sample, {
    req(nzchar(input$project_sample %||% ""))
    if (!file.exists(input$project_sample)) return(showNotification("The file of this sample is not found.", type = "error"))
    open_path(input$project_sample)
  })

  # samples not counted: those of the loaded file (all of them in a trio file), and the samples
  # of its index when the file is in the project
  exclude <- reactive({
    req(data())
    own <- sample_names(names(data()))
    file <- loaded()$file
    if (nzchar(input$project %||% "") && length(file)) {
      m <- index_info(input$project, file)
      own <- c(own, unlist(m$samples))
    }
    unique(own)
  })

  project <- reactive({
    id <- input$project
    if (!nzchar(id %||% "") || is.null(data())) return(NULL)
    indexed_files()
    p <- read_project(id)
    if (is.null(p)) return(NULL)
    counts <- tryCatch(project_counts(id, data(), exclude()), error = function(e) {
      showNotification(paste("Project:", conditionMessage(e)), type = "error"); NULL })
    if (is.null(counts)) return(NULL)
    list(id = id, name = p$name, counts = counts$counts, total = counts$total, exclude = exclude())
  })

  output$project_info <- renderUI({
    if (!nzchar(input$project %||% "")) return(helpText("Choose the project of this sample: each variant says in how many other samples of the project it is (columns In project / Hom in project)."))
    s <- project_state()
    req(s)
    st <- s$status
    pending <- sum(!st$indexed & st$exists)
    helpText(sprintf("%d files, %d samples indexed%s%s.", nrow(st), length(unique(unlist(strsplit(st$samples[st$indexed], ", ", fixed = TRUE)))),
                     if (pending) sprintf(", %d pending", pending) else "",
                     if (identical(s$indexing$stage, "running")) sprintf(" (indexing %d of %d...)", s$indexing$done + 1, s$indexing$total) else ""),
             if (!is.null(project())) sprintf("This sample is compared with %d other samples.", project()$total))
  })

  # ---- Manage projects (dialog) ----
  managed <- reactiveVal("")
  shinyFiles::shinyFileChoose(input, "project_add", roots = file_roots())
  observeEvent(input$manage_projects, {
    managed(if (nzchar(input$project %||% "")) input$project else (project_choices()[1] %||% ""))
    showModal(modalDialog(
      title = "Projects", size = "l", easyClose = TRUE, footer = modalButton("Close"),
      fluidRow(
        column(6, selectInput(ns("managed_project"), "Project", choices = project_choices(), selected = managed())),
        column(6, textInput(ns("new_project_name"), "New project", placeholder = "Name (e.g. Retina cohort 2026)"),
                  actionButton(ns("create_project"), "Create", icon = icon("plus")))),
      uiOutput(ns("managed_files")),
      div(style = "margin-top: 10px;",
          shinyFiles::shinyFilesButton(ns("project_add"), "Add variant files", "Choose the variant files of the project (several with Ctrl / Shift)",
                                       multiple = TRUE, icon = icon("folder-open")),
          actionButton(ns("index_project"), "Index the new files", icon = icon("sync")),
          actionLink(ns("delete_project"), "Delete the project", icon = icon("trash"), style = "margin-left: 15px; color: #c62828;")),
      helpText("One file per sample (a file with several samples counts each sample by its <sample>_GT column). Each file is",
               "read once in the background (big genomes take some minutes) and only its variants and genotypes are kept,",
               "in the shared folder; changed files are read again.")))
  })
  observeEvent(input$managed_project, managed(input$managed_project), ignoreInit = TRUE)
  managed_state <- reactive({ changed(); id <- managed(); req(nzchar(id)); state_of(id) })

  output$managed_files <- renderUI({
    if (!nzchar(managed())) return(helpText("There are no projects yet: create one."))
    s <- managed_state()
    st <- s$status
    ix <- s$indexing
    tagList(
      if (identical(ix$stage, "running")) div(class = "text-info", icon("spinner", class = "fa-spin"),
                                              sprintf(" Indexing %s (%d of %d)...", ix$file, ix$done + 1, ix$total)),
      if (length(unlist(ix$errors))) div(class = "text-danger", "Not indexed: ", paste(unlist(ix$errors), collapse = "; ")),
      if (!nrow(st)) helpText("No files yet: add the variant files of the samples.")
      else tags$table(class = "table table-condensed", style = "font-size: 13px;",
        tags$thead(tags$tr(tags$th("File"), tags$th("Samples"), tags$th("Index"), tags$th(""))),
        tags$tbody(lapply(seq_len(nrow(st)), function(i) tags$tr(
          tags$td(title = st$path[i], basename(st$path[i])),
          tags$td(st$samples[i]),
          tags$td(if (!st$exists[i]) span(class = "text-danger", "file not found")
                  else if (st$indexed[i]) span(class = "text-success", icon("check"), "indexed")
                  else span(class = "text-warning", "pending")),
          tags$td(tags$a(href = "#", title = "Remove from the project", icon("times"),
                         onclick = sprintf("Shiny.setInputValue('%s', %s, {priority: 'event'}); return false;",
                                           ns("remove_project_file"), jsonlite::toJSON(st$path[i], auto_unbox = TRUE)))))))))
  })

  observeEvent(input$create_project, {
    id <- tryCatch(create_project(input$new_project_name, user()$id), error = function(e) e)
    if (inherits(id, "error")) return(showNotification(conditionMessage(id), type = "error"))
    updateTextInput(session, "new_project_name", value = "")
    updateSelectInput(session, "managed_project", choices = project_choices(), selected = id)
    managed(id)
    refresh_choices()
  })
  observeEvent(input$project_add, {
    value <- input$project_add
    if (is.null(value) || is.numeric(value) || !nzchar(managed())) return()
    paths <- as.character(shinyFiles::parseFilePaths(file_roots(), value)$datapath)
    if (!length(paths)) return()
    n <- add_project_files(managed(), paths, user()$id)
    changed(changed() + 1)
    if (n) start_project_indexing(managed())
    showNotification(sprintf("%d files added%s.", n, if (n) "; indexing them in the background" else ""), type = "message")
  })
  observeEvent(input$index_project, {
    req(nzchar(managed()))
    if (!start_project_indexing(managed())) showNotification("The project is already being indexed.", type = "warning")
    changed(changed() + 1)
  })
  observeEvent(input$remove_project_file, {
    req(nzchar(managed()))
    remove_project_file(managed(), input$remove_project_file)
    changed(changed() + 1)
  })
  observeEvent(input$delete_project, {
    req(nzchar(managed()))
    p <- read_project(managed())
    if (!identical(p$created_by, user()$id) && !isTRUE(user()$is_admin))
      return(showNotification("Only who created the project or an administrator can delete it.", type = "warning"))
    shinyalert::shinyalert("Delete the project?", sprintf("%s and its index (the variant files are not touched).", p$name),
                           type = "warning", showCancelButton = TRUE, confirmButtonText = "Delete",
                           callbackR = function(ok) if (isTRUE(ok)) {
                             delete_project(managed())
                             managed(project_choices()[1] %||% "")
                             updateSelectInput(session, "managed_project", choices = project_choices(), selected = managed())
                             refresh_choices(selected = "")
                           })
  })

  project
}
