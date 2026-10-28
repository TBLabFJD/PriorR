# PanelApp of the SNV tab: each user can download the latest version for themselves
# (shared/panelapp/<user>, panels.R); the panel lists follow the copy of the user.

snv_panelapp_server <- function(input, output, session, user) {
  ns <- session$ns
  panelapp_job <- reactiveValues(running = FALSE, log = NULL, done = 0)
  refresh_panel_choices <- function() {
    updateSelectizeInput(session, "panels", choices = panel_choices(specific = FALSE), selected = isolate(input$panels), server = TRUE)
    updateSelectizeInput(session, "specific_panels", choices = panel_choices(specific = TRUE), selected = isolate(input$specific_panels), server = TRUE)
    updateSelectizeInput(session, "excluded_panels", choices = panel_choices(), selected = isolate(input$excluded_panels), server = TRUE)
    updateSelectizeInput(session, "glow_run_panel", choices = panel_choices(min_genes = glow_min_genes), selected = isolate(input$glow_run_panel))
  }
  observeEvent(user(), refresh_panel_choices())
  output$red_note <- renderUI({
    req(isTRUE(input$red), user())
    panelapp_job$done
    if (!panelapp_has_red())
      div(style = "margin: 0 15px 6px 15px; font-size: 12px; color: #ffcc80; white-space: normal;", icon("exclamation-triangle"),
          "This PanelApp copy has no red genes (made before PriorR 4.0): download the latest PanelApp below.")
  })
  output$panelapp_info <- renderUI({
    req(user())
    panelapp_job$done
    info <- panelapp_info()
    div(style = "margin: 10px 15px 4px 15px; font-size: 12px; white-space: normal;",
        if (panelapp_job$running) tagList(icon("spinner", class = "fa-spin"), "Downloading PanelApp (10-15 minutes)...")
        else sprintf("PanelApp of %s (%s copy).", info$date, if (info$own) "your own" else "the lab"))
  })
  observeEvent(input$update_panelapp, {
    req(user(), !panelapp_job$running)
    panelapp_job$log <- tempfile(fileext = ".log")
    system2(file.path(R.home("bin"), "Rscript"),
            c(shQuote(normalizePath(file.path(data_dir, "scripts", "update_panels.R"))), shQuote(user_panelapp_dir(user()$username))),
            stdout = panelapp_job$log, stderr = panelapp_job$log, wait = FALSE)
    panelapp_job$running <- TRUE
  })
  observe({
    req(panelapp_job$running)
    invalidateLater(5000)
    log <- if (file.exists(panelapp_job$log)) readLines(panelapp_job$log, warn = FALSE) else character(0)
    if (any(grepl("^Done", log))) {
      panelapp_job$running <- FALSE
      panelapp_job$done <- panelapp_job$done + 1
      refresh_panel_choices()
      showNotification("PanelApp updated (only for you)", type = "message")
    } else if (any(grepl("Error|Execution halted", log))) {
      panelapp_job$running <- FALSE
      showNotification(paste(c("PanelApp could not be updated:", tail(log, 2)), collapse = " "), type = "error", duration = NULL)
    }
  })
  invisible(refresh_panel_choices)
}
