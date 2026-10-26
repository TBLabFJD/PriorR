# First-run setup on a computer: where is the shared folder of the lab? It's asked only when
# PriorR runs on the user's computer and no priorr.cfg (of the user or of the lab) says it.

setup_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    roots <- c(Home = path.expand("~"), shinyFiles::getVolumes()())
    shinyFiles::shinyDirChoose(input, "folder", roots = roots)
    chosen <- reactive({
      if (is.null(input$folder) || is.numeric(input$folder)) return(NULL)
      path <- shinyFiles::parseDirPath(roots, input$folder)
      if (length(path)) as.character(path) else NULL
    })

    showModal(modalDialog(
      title = paste("Set up", priorr_name, "on this computer"), easyClose = FALSE, footer = NULL, size = "m",
      div(style = "text-align: center; margin-bottom: 15px;", tags$img(src = "logo_v4.png", style = "width: 200px;")),
      p("PriorR keeps the users, comments, profiles and saved sessions in a ",
        strong("shared folder of the lab"), ", so that every computer sees the same. ",
        "The variant files are not copied: they are opened from where they are."),
      radioButtons(ns("mode"), NULL, choices = c(
        "Use a shared folder of the lab" = "shared",
        "Use PriorR only on this computer" = "local")),
      conditionalPanel(sprintf("input['%s'] == 'shared'", ns("mode")),
        shinyFiles::shinyDirButton(ns("folder"), "Choose the shared folder", "Shared folder of the lab",
                                   icon = icon("folder-open")),
        div(style = "margin: 6px 0;", textOutput(ns("folder_path"))),
        helpText("If you are the first person of the lab, choose an empty folder of the shared drive;",
                 "everybody else must choose the same one.")),
      conditionalPanel(sprintf("input['%s'] == 'local'", ns("mode")),
        helpText("The data will be kept in:", local_shared_dir(),
                 "It can be changed later in the configuration file (see Help).")),
      uiOutput(ns("error")),
      actionButton(ns("save"), "Save and continue", class = "btn-success", style = "width: 100%;")))

    output$folder_path <- renderText(if (is.null(chosen())) "No folder chosen" else chosen())

    observeEvent(input$save, {
      folder <- if (input$mode == "local") local_shared_dir() else chosen()
      if (is.null(folder)) return(output$error <- renderUI(p(class = "text-danger", "Choose the shared folder")))
      if (input$mode == "shared" && !dir.exists(folder))
        return(output$error <- renderUI(p(class = "text-danger", "The folder doesn't exist")))
      dir.create(folder, showWarnings = FALSE, recursive = TRUE)
      if (file.access(folder, 2) != 0)
        return(output$error <- renderUI(p(class = "text-danger", "PriorR can't write in that folder")))
      save_user_config(list(shared_folder = folder))
      session$reload()
    })
  })
}

# Data folder when PriorR is used only on this computer
local_shared_dir <- function() Sys.getenv("PRIORR_LOCAL_SHARED", file.path(tools::R_user_dir("priorr", "data"), "shared"))
