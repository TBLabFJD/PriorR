# Help page (help/help.md) with the configuration of this computer.

tblab_url <- "https://www.translationalbioinformaticslab.es/"
# Who develops PriorR
priorr_team <- "Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación Jiménez Díaz"

help_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(8, div(style = "max-width: 900px; font-size: 14px;",
                  HTML(commonmark::markdown_html(paste(readLines("help/help.md", encoding = "UTF-8"), collapse = "\n"),
                                                 extensions = TRUE)))),
    column(4,
      shinydashboardPlus::box(
        title = "About PriorR", width = NULL, status = "primary", collapsible = FALSE,
        div(style = "display: flex; align-items: center; gap: 12px;",
            tags$a(href = tblab_url, target = "_blank", tags$img(src = "tblab_logo.png", style = "width: 70px;", alt = "Bioinformatics Unit")),
            div(strong(priorr_name), br(),
                "Bioinformatics Unit", br(),
                "Genetics and Genomics Department", br(),
                "IIS-Fundación Jiménez Díaz")),
        tags$ul(style = "margin-top: 10px; padding-left: 18px;",
          tags$li(tags$a(href = tblab_url, target = "_blank", "Website of the Bioinformatics Unit")),
          tags$li(tags$a(href = "https://github.com/TBLabFJD/PriorR", target = "_blank", "PriorR on GitHub")),
          tags$li(tags$a(href = "https://github.com/TBLabFJD/GLOWgenes", target = "_blank", "GLOWgenes on GitHub")),
          tags$li(tags$a(href = "PriorR4_UserGuide.pdf", target = "_blank", strong(paste("PriorR", priorr_version, "user guide (PDF)"))))),
        helpText("GLOWgenes: de la Fuente L et al., Int J Mol Sci 2023;24(2):1661.", br(),
                 "GLOWgenes rankings of PanelApp, disease-specific genes and SGDS: Benítez Y, Uría-Regojo G, Mínguez P. A global survey of systems biology-based predictions of gene-rare disease associations to enhance new diagnoses. Sci Rep 2026. doi:10.1038/s41598-026-54510-3 (CC BY-NC-ND 4.0).", br(),
                 "RetiGene: Am J Hum Genet 2025;112(10):2253-2265.")),
      shinydashboardPlus::box(
        title = "Configuration of this computer", width = NULL, status = "success", collapsible = FALSE,
        uiOutput(ns("config")))))
}

help_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$config <- renderUI({
      file_info <- function(path) if (file.exists(path)) path else paste(path, "(doesn't exist)")
      tags$dl(
        tags$dt("Shared folder (shared_folder)"), tags$dd(code(shared_dir)),
        tags$dt("Configuration of this computer"), tags$dd(code(file_info(user_config_file()))),
        tags$dt("Configuration of the lab"), tags$dd(code(file_info(lab_config_file()))),
        tags$dt("Cache of big files (cache_folder)"), tags$dd(code(wgs_cache_dir)),
        tags$dt("Genome build (genome_build)"), tags$dd(genome_build()),
        tags$dt("Memory / processors"), tags$dd(paste(wgs_memory_limit, "/", wgs_threads)),
        tags$dt("Version"), tags$dd(paste(priorr_name, "· R", getRversion())))
    })
  })
}
