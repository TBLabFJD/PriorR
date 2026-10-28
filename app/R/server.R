server <- function(input, output, session) {
  # First time on this computer and no configuration: where is the shared folder?
  if (!shared_folder_configured()) return(setup_server("setup"))

  user <- auth_server("auth")
  # the user of the session (their own PanelApp copy, panels.R); before the modules react to user()
  observeEvent(user(), session$userData$username <- user()$username, priority = 100)

  snv_server("snv", user)
  cnv_server("cnv", user)
  admin_server("admin", user)
  bookmarks_server("bookmarks", user)
  help_server("help")
  observeEvent(input$open_help, updateTabsetPanel(session, "navbar", selected = "help"))

  # The sidebar follows the selected tab; the administration tab only for administrators
  observeEvent(input$navbar, {
    shinyjs::toggle("snv_sidebar", condition = input$navbar != "cnv")
    shinyjs::toggle("cnv_sidebar", condition = input$navbar == "cnv")
  })
  hideTab("navbar", "admin")
  observeEvent(user(), if (isTRUE(user()$is_admin)) showTab("navbar", "admin"))

  setBookmarkExclude(c("navbar", "open_help"))
}
