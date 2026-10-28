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
  if (!interactive() && nzchar(Sys.getenv("PRIORR_ADMIN_PASSWORD"))) return(Sys.getenv("PRIORR_ADMIN_PASSWORD"))   # tests
  cat("Password: ")
  system("stty -echo", ignore.stderr = TRUE)
  on.exit(system("stty echo", ignore.stderr = TRUE))
  password <- readLines(file("stdin"), n = 1)
  cat("\n")
  password
}

# Users are in accounts.json of the shared folder (store.R); their id is the user name
user_id <- function(username) {
  a <- read_accounts()
  i <- find_user(a, username)
  if (is.na(i)) stop("No user ", username, " in ", shared_dir)
  a$users[[i]]$username
}

read_new_password <- function() {
  p1 <- read_password()
  cat("Again: ")
  p2 <- read_password()
  if (!identical(p1, p2)) stop("The passwords are not the same")
  p1
}

switch(action,
  "create-admin" = {
    create_user(username, read_new_password(), full_name = username, is_admin = TRUE)
    cat("Administrator", username, "created in", shared_dir, "\n")
  },
  "reset-password" = {
    id <- user_id(username)            # the user first, then the password
    set_password(id, read_new_password())
    cat("Password of", username, "changed (in", shared_dir, ")\n")
  },
  "import-v3" = {
    path <- if (length(args) >= 3) args[3] else file.path(data_dir, "Database", "variants_priorr.txt")
    n <- import_legacy_classifications(path, user_id(username))
    cat(n, "classifications imported as comments of", username, "\n")
  },
  stop("Unknown action ", action))
