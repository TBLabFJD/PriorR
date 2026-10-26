# "Open file": PriorR runs on the user's own computer, so big files are opened directly from
# the disk (local or shared folder) instead of being uploaded: no size limit and no copy.

# `sidebar`: shown in the dark sidebar (lighter text) or in a dialog
open_file_ui <- function(id, label, title = label, sidebar = TRUE) {
  tagList(
    shinyFiles::shinyFilesButton(id, label, title, multiple = FALSE, icon = icon("folder-open"),
                                 style = if (sidebar) "margin: 6px 15px 0 15px;" else ""),
    div(style = sprintf("margin: 4px %s 6px %s; font-size: 13px; color: %s; word-break: break-all;",
                        if (sidebar) "15px" else "0", if (sidebar) "15px" else "0", if (sidebar) "#b8c7ce" else "#333"),
        textOutput(paste0(id, "_path"), inline = TRUE)))
}

# Inside a module server: returns a reactive with the chosen path (NULL until a file is chosen).
# Every file of the folder is shown (annotation pipelines name their files in many ways); a file
# PriorR can't read gives an error when it is loaded. `filetypes` only limits the list when given.
open_file_server <- function(id, input, output, filetypes = NULL) {
  roots <- file_roots()
  if (is.null(filetypes)) shinyFiles::shinyFileChoose(input, id, roots = roots)
  else shinyFiles::shinyFileChoose(input, id, roots = roots, filetypes = filetypes)
  path <- reactive({
    value <- input[[id]]
    if (is.null(value) || is.numeric(value)) return(NULL)
    chosen <- shinyFiles::parseFilePaths(roots, value)
    if (nrow(chosen) == 0) NULL else as.character(chosen$datapath[1])
  })
  output[[paste0(id, "_path")]] <- renderText(if (is.null(path())) "" else basename(path()))
  path
}
