# Downloads the PanelApp panels (Genomics England and Australia) into Program_files/Panels.
# Run from the app/ folder:  Rscript Program_files/scripts/update_panels.R [folder]
# Without folder: the lab copy (Program_files/Panels); PriorR passes the folder of a user.
suppressMessages(library(shiny))
for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else panels_dir
update_panelapp(function(source, i, n) if (i %% 50 == 0 || i == n) message(source, ": ", i, " / ", n), dir = dir)
message("Done: ", paste(list.files(dir), collapse = ", "))
