# Administration from the command line (run from the app/ folder):
#
#   Rscript Program_files/scripts/priorr_admin.R create-admin <user>      # asks for the password
#   Rscript Program_files/scripts/priorr_admin.R reset-password <user>    # asks for the password
#   Rscript Program_files/scripts/priorr_admin.R import-v3 <user> [file]  # PriorR 3 classifications
#
# Users and comments are in the shared folder (see store.R); import-v3 reads Database/variants_priorr.txt of PriorR 3.

suppressMessages(library(shiny))
for (f in list.files("R", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: priorr_admin.R create-admin|reset-password|import-v3 <user> [file]")
action <- args[1]
username <- args[2]

read_password <- function() {
  cat("Password: ")
  system("stty -echo", ignore.stderr = TRUE)
  on.exit(system("stty echo", ignore.stderr = TRUE))
  password <- readLines(file("stdin"), n = 1)
  cat("\n")
  password
}

user_id <- function(username) {
  id <- db_query("SELECT id FROM users WHERE username = ?", username)$id
  if (length(id) != 1) stop("No user ", username)
  id
}

switch(action,
  "create-admin" = {
    create_user(username, read_password(), is_admin = TRUE)
    cat("Administrator", username, "created in", db_path, "\n")
  },
  "reset-password" = {
    set_password(user_id(username), read_password())
    cat("Password of", username, "changed\n")
  },
  "import-v3" = {
    path <- if (length(args) >= 3) args[3] else file.path(data_dir, "Database", "variants_priorr.txt")
    n <- import_legacy_classifications(path, user_id(username))
    cat(n, "classifications imported as comments of", username, "\n")
  },
  stop("Unknown action ", action))
