# Page layout.

ui <- function(request) {
  shinydashboardPlus::dashboardPage(
    skin = "green",
    title = priorr_name,

    header = shinydashboardPlus::dashboardHeader(
      title = tags$img(src = "logo_header_v4.png", height = "40px", alt = priorr_name),
      titleWidth = 280,
      leftUi = tagList(
        tags$li(class = "dropdown",
                actionLink("open_help", "Help", icon = icon("question-circle"), style = "color: white;")),
        tags$li(class = "dropdown", snv_patient_ui("snv"))),
      auth_header_ui("auth")),

    sidebar = shinydashboardPlus::dashboardSidebar(
      width = 280,
      div(id = "snv_sidebar", snv_sidebar_ui("snv")),
      shinyjs::hidden(div(id = "cnv_sidebar", cnv_sidebar_ui("cnv"))),
      tags$a(href = tblab_url, target = "_blank", title = priorr_team,
             style = "display: flex; align-items: center; gap: 8px; margin: 25px 15px 15px 15px; color: #b8c7ce; font-size: 12px;",
             tags$img(src = "tblab_logo.png", style = "width: 34px; border-radius: 4px;"), "PriorR \u00b7 Bioinformatics Unit, IIS-FJD")),

    body = shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),
      table_css,
      tags$head(
        tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
        tags$style(HTML("
          .main-header .logo, .main-header .logo:hover { padding: 0 !important; }
          .sidebar-menu .help-block, .sidebar-menu .shiny-html-output { white-space: normal; margin-left: 15px; margin-right: 15px; }
          /* long labels of the sidebar wrap instead of going into the main panel */
          .main-sidebar .pretty { white-space: normal; width: calc(100% - 30px); }
          .main-sidebar .pretty .state label { text-indent: 0; padding-left: 1.7em; line-height: 1.3; }
          .main-sidebar .pretty .state label:before, .main-sidebar .pretty .state label:after { top: 0.05em; }
          .main-sidebar .action-button, .main-sidebar .shiny-download-link { white-space: normal; max-width: calc(100% - 30px); text-align: left; }
          .main-sidebar .sidebar-menu > li > a > span { display: inline-block; max-width: 200px; white-space: normal; vertical-align: middle; }
          .main-sidebar .irs { margin-right: 12px; }
          .main-sidebar .radio label, .main-sidebar .checkbox label, .main-sidebar .control-label { white-space: normal; }
          /* view buttons above the table: the view not in use in grey */
          /* lists of the sidebar (columns...): the list scrolls inside a menu that never cuts it */
          .main-sidebar .bootstrap-select .dropdown-menu.open { max-height: none !important; overflow: visible !important; }
          .main-sidebar .bootstrap-select div.inner { max-height: 300px !important; overflow-y: auto !important; }
          .main-sidebar .bootstrap-select ul.dropdown-menu.inner { max-height: none !important; overflow: visible !important; padding-bottom: 4px; }
          .pr-viewbar > :first-child { flex-shrink: 0; }
          .pr-viewbar #snv-view .btn-group { display: inline-flex; }
          #snv-view .btn:not(.active) { background: #e0e0e0; border-color: #cfcfcf; color: #757575; box-shadow: none; }
          #snv-view .btn:not(.active):hover { background: #d5d5d5; color: #424242; }
          .pr-patient { display: inline-block; margin: 9px 0 0 25px; padding: 5px 14px; border-radius: 16px;
                        background: rgba(255, 255, 255, 0.95); color: #1b5e20; font-size: 16px; white-space: nowrap; }
          .pr-patient b { font-size: 18px; margin-left: 4px; letter-spacing: 0.3px; }
          .pr-patient .fa { margin-right: 6px; }
          .pr-patient-others { color: #607d8b; font-size: 13px; margin-left: 8px; }
          .main-header .logo img { display: block; height: 40px; margin: 5px auto; }"))),
      tabsetPanel(
        id = "navbar",
        tabPanel("SNV", value = "snv", snv_main_ui("snv")),
        tabPanel("CNV", value = "cnv", cnv_main_ui("cnv")),
        tabPanel("Users and groups", value = "admin", br(), admin_ui("admin")),
        tabPanel("Help", value = "help", icon = icon("question-circle"), br(), help_ui("help"))),
      br(),
      bookmarks_ui("bookmarks")))
}
