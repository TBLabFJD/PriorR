# Starts PriorR on this computer and opens it in the browser.
# Started by PriorR.bat (Windows), PriorR.command (macOS) or priorr.sh (Linux).
#
# The first time on each operating system it installs the R packages that PriorR needs in
# library/<platform>/ next to this file (in the shared folder), so the next computers with the
# same system don't install anything.

args <- commandArgs(FALSE)
priorr_home <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))

packages <- c("shiny", "DT", "dplyr", "data.table", "jsonlite", "ggplot2", "shinydashboard",
              "shinydashboardPlus", "shinyWidgets", "shinyalert", "shinyjs", "shinyFiles", "sodium",
              "DBI", "duckdb", "commonmark")

platform <- paste0(R.version$platform, "-R", getRversion()$major, ".", getRversion()$minor)
lib <- file.path(priorr_home, "library", platform)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
# PriorR installed where the user can't write (e.g. by an administrator): the packages go to the
# folder of the user
if (file.access(lib, 2) != 0) {
  lib <- file.path(tools::R_user_dir("priorr", "data"), "library", platform)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(lib, .libPaths()))
# The R processes that PriorR starts in the background (loading big files, PanelApp, GLOWgenes) must
# find the same packages: they inherit R_LIBS
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("PriorR: installing ", length(missing), " R packages (only the first time on this system)...")
  repos <- "https://cloud.r-project.org"
  # Linux: binary packages from Posit Package Manager when the distribution is known
  if (Sys.info()[["sysname"]] == "Linux" && file.exists("/etc/os-release")) {
    os <- readLines("/etc/os-release")
    codename <- sub("^VERSION_CODENAME=", "", grep("^VERSION_CODENAME=", os, value = TRUE))
    if (length(codename) && nzchar(codename)) {
      repos <- c(sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename), repos)
      options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                      paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))
    }
  }
  install.packages(missing, lib = lib, repos = repos)
}

# Running on the user's own computer: PriorR can remember the user ("Remember me")
if (Sys.getenv("PRIORR_LOCAL") == "") Sys.setenv(PRIORR_LOCAL = "1")

# Always the same port (8484, or PRIORR_PORT), so the browser keeps "Remember me" and the address
# doesn't change; when it is in use (e.g. another PriorR), the next free one
port_free <- function(p) tryCatch({
  s <- httpuv::startServer("127.0.0.1", p, list())
  httpuv::stopServer(s)
  TRUE
}, error = function(e) FALSE)
port <- as.integer(Sys.getenv("PRIORR_PORT", "8484"))
if (!port_free(port)) {
  free <- Filter(port_free, port + 1:20)
  port <- if (length(free)) free[1] else httpuv::randomPort()
}
message("PriorR v4.0: http://127.0.0.1:", port, "  (close this window to stop PriorR)")
shiny::runApp(file.path(priorr_home, "app"), host = "127.0.0.1", port = port,
              launch.browser = Sys.getenv("PRIORR_NO_BROWSER") == "")
