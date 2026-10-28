# PriorR: prioritization of disease-linked genetic variants.
# Run with shiny::runApp("app"); the files in R/ are loaded automatically.

library(shiny)

# Saved sessions are stored in the shared folder (with their uploaded files), so they can be
# opened from any computer
bookmark_dir <- function() Sys.getenv("PRIORR_BOOKMARK_DIR", file.path(shared_dir, "bookmarks"))
shinyOptions(
  save.interface = function(id, callback) {
    path <- file.path(bookmark_dir(), id)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    callback(path)
  },
  load.interface = function(id, callback) callback(file.path(bookmark_dir(), id)))

enableBookmarking(store = "server")

shinyApp(ui, server)
