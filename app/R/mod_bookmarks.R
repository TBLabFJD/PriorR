# Saved sessions (Shiny bookmarks) of each user.

bookmarks_ui <- function(id) {
  ns <- NS(id)
  shinydashboardPlus::box(
    title = "Saved sessions", width = NULL, status = "warning", collapsible = TRUE, collapsed = TRUE,
    fluidRow(
      column(6, textInput(ns("description"), NULL, placeholder = "Description of the analysis")),
      column(4, actionButton(ns("save"), "Save session", icon = icon("bookmark"), class = "btn-warning"))),
    DT::DTOutput(ns("list")))
}

bookmarks_server <- function(id, user) {
  moduleServer(id, function(input, output, session) {
    changed <- reactiveVal(0)
    pending <- reactiveVal(NULL)

    # onBookmarked() isn't available inside modules: use the root session
    root <- session$rootScope()

    observeEvent(input$save, {
      req(user())
      pending(input$description)
      root$doBookmark()
    })

    root$onBookmarked(function(url) {
      req(user(), !is.null(pending()))
      add_bookmark(user()$id, pending(), sub(".*[?&]_state_id_=([^&]+).*", "\\1", url))
      pending(NULL)
      updateTextInput(session, "description", value = "")
      changed(changed() + 1)
      showNotification("Session saved", type = "message")
    })

    output$list <- DT::renderDT({
      req(user())
      changed()
      b <- list_bookmarks(user()$id)
      b$state_id <- sprintf("<a href='?_state_id_=%s'>Open</a>", b$state_id)
      DT::datatable(b, escape = FALSE, rownames = FALSE, selection = "none",
                    colnames = c("Description", "Session", "Saved"), options = list(dom = "tp", pageLength = 5))
    })

    # The inputs of this box are not part of the saved analysis
    setBookmarkExclude(c("description", "save", "list_rows_all", "list_rows_current", "list_rows_selected",
                         "list_search", "list_state", "list_cell_clicked"))
  })
}
