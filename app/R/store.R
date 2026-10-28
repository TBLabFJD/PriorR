# Shared data of PriorR (users, groups, comments, saved sessions) as plain files in a shared
# folder that every computer running PriorR mounts (PRIORR_SHARED_DIR):
#
#   accounts.json            users and groups (written with a lock: administrators, password changes)
#   comments/<user>.jsonl    comments of each user (only that user writes the file)
#   sessions/<user>.jsonl    saved sessions of each user (only that user writes the file)
#   bookmarks/<id>/          state of the saved sessions
#
# No file is ever written by two computers at the same time, so it's safe on network folders.
#
# Comment visibility:
#   - "private": only the author
#   - "groups":  the author and every user sharing a group with the author
#   - "group":   the author and the members of one group

store_file <- function(...) file.path(shared_dir, ...)

user_file_name <- function(username) paste0(gsub("[^A-Za-z0-9._-]", "_", tolower(username)), ".jsonl")

now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# ---- Accounts (users and groups) ----

empty_accounts <- function() list(users = list(), groups = list())

read_accounts <- function() {
  f <- store_file("accounts.json")
  if (!file.exists(f)) return(empty_accounts())
  jsonlite::fromJSON(f, simplifyVector = FALSE)
}

# Read-modify-write of accounts.json under a lock (a folder: creating it is atomic)
update_accounts <- function(change) {
  dir.create(shared_dir, showWarnings = FALSE, recursive = TRUE)
  lock <- store_file("accounts.lock")
  for (i in 1:100) {
    if (dir.create(lock, showWarnings = FALSE)) break
    info <- file.info(lock)
    if (!is.na(info$mtime) && difftime(Sys.time(), info$mtime, units = "secs") > 60) unlink(lock, recursive = TRUE)
    Sys.sleep(0.1)
    if (i == 100) stop("The accounts file is locked, try again")
  }
  on.exit(unlink(lock, recursive = TRUE))
  accounts <- change(read_accounts())
  tmp <- store_file(sprintf("accounts.json.%s.tmp", Sys.getpid()))
  jsonlite::write_json(accounts, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null")
  file.rename(tmp, store_file("accounts.json"))
  invisible(accounts)
}

find_user <- function(accounts, username) {
  i <- which(vapply(accounts$users, function(u) tolower(u$username) == tolower(trimws(username)), logical(1)))
  if (length(i)) i[1] else NA_integer_
}

user_count <- function() length(read_accounts()$users)

check_password <- function(password) {
  if (nchar(password) < 8) stop("The password must have at least 8 characters")
}

# Returns the user id (its user name)
create_user <- function(username, password, full_name = "", is_admin = FALSE) {
  username <- trimws(username)
  if (!nzchar(username)) stop("The user name can't be empty")
  check_password(password)
  update_accounts(function(a) {
    if (!is.na(find_user(a, username))) stop("The user ", username, " already exists")
    a$users[[length(a$users) + 1]] <- list(
      username = username, full_name = full_name, password_hash = sodium::password_store(password),
      is_admin = isTRUE(is_admin), active = TRUE, groups = list(), created_at = now())
    a
  })
  username
}

# The user (id, username, full_name, is_admin) or NULL
authenticate <- function(username, password) {
  a <- read_accounts()
  i <- find_user(a, username)
  if (is.na(i)) return(NULL)
  u <- a$users[[i]]
  if (!isTRUE(u$active) || !sodium::password_verify(u$password_hash, password)) return(NULL)
  list(id = u$username, username = u$username, full_name = u$full_name, is_admin = isTRUE(u$is_admin))
}

modify_user <- function(user_id, change) {
  update_accounts(function(a) {
    i <- find_user(a, user_id)
    if (is.na(i)) stop("No user ", user_id)
    a$users[[i]] <- change(a$users[[i]])
    a
  })
}

set_password <- function(user_id, password) {
  check_password(password)
  modify_user(user_id, function(u) { u$password_hash <- sodium::password_store(password); u$remember <- list(); u })
}

update_user <- function(user_id, full_name = NULL, is_admin = NULL, active = NULL) {
  modify_user(user_id, function(u) {
    if (!is.null(full_name)) u$full_name <- full_name
    if (!is.null(is_admin)) u$is_admin <- isTRUE(is_admin)
    if (!is.null(active)) u$active <- isTRUE(active)
    u
  })
}

list_users <- function() {
  users <- read_accounts()$users
  if (!length(users)) return(data.frame(id = character(0), username = character(0), full_name = character(0),
                                        is_admin = logical(0), active = logical(0), created_at = character(0), groups = character(0)))
  data.frame(
    id = vapply(users, `[[`, "", "username"),
    username = vapply(users, `[[`, "", "username"),
    full_name = vapply(users, function(u) u$full_name %||% "", ""),
    is_admin = vapply(users, function(u) isTRUE(u$is_admin), logical(1)),
    active = vapply(users, function(u) isTRUE(u$active), logical(1)),
    created_at = vapply(users, function(u) u$created_at %||% "", ""),
    groups = vapply(users, function(u) paste(sort(unlist(u$groups)), collapse = ", "), ""),
    stringsAsFactors = FALSE)
}

# ---- Groups (identified by their name) ----

create_group <- function(name, description = "") {
  name <- trimws(name)
  if (!nzchar(name)) stop("The group name can't be empty")
  update_accounts(function(a) {
    if (tolower(name) %in% tolower(vapply(a$groups, `[[`, "", "name"))) stop("The group ", name, " already exists")
    a$groups[[length(a$groups) + 1]] <- list(name = name, description = description)
    a
  })
  name
}

list_groups <- function() {
  a <- read_accounts()
  names <- vapply(a$groups, `[[`, "", "name")
  members <- vapply(names, function(g) sum(vapply(a$users, function(u) g %in% unlist(u$groups), logical(1))), 0L)
  data.frame(id = names, name = names, description = vapply(a$groups, function(g) g$description %||% "", ""),
             members = unname(members), stringsAsFactors = FALSE)
}

groups_of <- function(accounts, user_id) {
  i <- find_user(accounts, user_id)
  if (is.na(i)) character(0) else as.character(unlist(accounts$users[[i]]$groups))
}

user_groups <- function(user_id) {
  g <- sort(groups_of(read_accounts(), user_id))
  data.frame(id = g, name = g, stringsAsFactors = FALSE)
}

set_user_groups <- function(user_id, group_ids) {
  modify_user(user_id, function(u) { u$groups <- as.list(unique(as.character(group_ids))); u })
}

# ---- Per-user append-only files ----

append_record <- function(folder, user_id, record) {
  dir <- store_file(folder)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  cat(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null"), "\n", sep = "",
      file = file.path(dir, user_file_name(user_id)), append = TRUE)
}

.files_cache <- new.env()

# Records of one file (cached while the file doesn't change; incomplete lines are skipped)
read_records <- function(path) {
  info <- file.info(path)
  key <- paste(path, info$size, info$mtime)
  if (!is.null(.files_cache[[path]]) && identical(.files_cache[[path]]$key, key)) return(.files_cache[[path]]$records)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  records <- lapply(lines[nzchar(lines)], function(l) tryCatch(jsonlite::fromJSON(l, simplifyVector = FALSE), error = function(e) NULL))
  records <- Filter(Negate(is.null), records)
  .files_cache[[path]] <- list(key = key, records = records)
  records
}

# ---- Comments ----

add_comment <- function(user_id, variant_key, gene = NA, classification = NA, comment = "",
                        visibility = "groups", group_id = NA, acmg = NULL, snapshot = NULL) {
  if (identical(classification, "")) classification <- NA
  if (is.na(classification) && !nzchar(trimws(comment %||% ""))) stop("Write a comment or choose a classification")
  if (!visibility %in% c("private", "groups", "group")) stop("Unknown visibility ", visibility)
  if (visibility == "group" && (is.null(group_id) || is.na(group_id))) stop("Choose the group to share the comment with")
  id <- paste0(user_file_name(user_id), "-", format(Sys.time(), "%Y%m%d%H%M%OS3"), "-", paste(sample(c(0:9, letters[1:6]), 6, TRUE), collapse = ""))
  append_record("comments", user_id, list(
    id = id, variant_key = variant_key, gene = if (is.na(gene)) NULL else gene,
    classification = if (is.na(classification)) NULL else classification, comment = comment,
    author = user_id, visibility = visibility, group = if (visibility == "group") as.character(group_id) else NULL,
    acmg = if (length(acmg)) as.list(acmg) else NULL,
    # the annotation when it was classified (reanalysis, acmg.R)
    snapshot = if (!is.na(classification) && length(snapshot)) snapshot else NULL,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")))
  invisible(id)
}

# Only the author can delete a comment (a deletion record in the author's file)
delete_comment <- function(comment_id, user_id) {
  mine <- all_comments()
  if (!comment_id %in% mine$id[tolower(mine$author) == tolower(user_id)]) return(FALSE)
  append_record("comments", user_id, list(deleted = comment_id))
  TRUE
}

all_comments <- function() {
  files <- list.files(store_file("comments"), pattern = "\\.jsonl$", full.names = TRUE)
  records <- unlist(lapply(files, read_records), recursive = FALSE)
  deleted <- unlist(lapply(records, `[[`, "deleted"))
  records <- Filter(function(r) is.null(r$deleted) && !(r$id %in% deleted), records)
  field <- function(name) vapply(records, function(r) if (is.null(r[[name]])) NA_character_ else as.character(r[[name]]), "")
  out <- data.frame(id = field("id"), variant_key = field("variant_key"), gene = field("gene"),
                    classification = field("classification"), comment = field("comment"), created_at = field("created_at"),
                    visibility = field("visibility"), group_name = field("group"), author = field("author"),
                    acmg = vapply(records, function(r) paste(unlist(r$acmg), collapse = ", "), ""),
                    stringsAsFactors = FALSE)
  out$snapshot <- lapply(records, function(r) r$snapshot)
  out
}

# Comments the user can see (of one variant, or all), newest first
visible_comments <- function(user_id, variant_key = NULL) {
  cm <- all_comments()
  if (!is.null(variant_key)) cm <- cm[cm$variant_key == variant_key, ]
  accounts <- read_accounts()
  my_groups <- groups_of(accounts, user_id)
  shares_group <- vapply(cm$author, function(a) length(intersect(groups_of(accounts, a), my_groups)) > 0, logical(1))
  visible <- tolower(cm$author) == tolower(user_id) |
    (cm$visibility == "groups" & shares_group) |
    (cm$visibility == "group" & cm$group_name %in% my_groups)
  cm <- cm[visible, ]
  cm$user_id <- cm$author
  cm[order(cm$created_at, decreasing = TRUE), ]
}

# Per variant: latest visible classification and number of visible comments
comment_annotations <- function(user_id) {
  cm <- visible_comments(user_id)
  if (nrow(cm) == 0) return(data.frame(variant_key = character(0), classification = character(0), n_comments = integer(0)))
  n <- table(cm$variant_key)
  classified <- cm[!is.na(cm$classification), ]          # newest first
  latest <- classified[!duplicated(classified$variant_key), ]
  data.frame(variant_key = names(n), n_comments = as.integer(n),
             classification = latest$classification[match(names(n), latest$variant_key)],
             stringsAsFactors = FALSE)
}

# ---- Saved sessions ----

# `state_id` is the Shiny bookmark id; the link is relative, so it works whatever the port is
add_bookmark <- function(user_id, description, state_id) {
  append_record("sessions", user_id, list(description = description, state_id = state_id, created_at = now()))
}

list_bookmarks <- function(user_id) {
  f <- file.path(store_file("sessions"), user_file_name(user_id))
  records <- if (file.exists(f)) read_records(f) else list()
  df <- data.frame(description = vapply(records, function(r) r$description %||% "", ""),
                   state_id = vapply(records, function(r) r$state_id %||% "", ""),
                   created_at = vapply(records, function(r) r$created_at %||% "", ""), stringsAsFactors = FALSE)
  df[rev(seq_len(nrow(df))), ]
}

# ---- Classifications of PriorR 3.x (Database/variants_priorr.txt) ----

# Imported as comments of `user_id` visible to its groups, marked as coming from PriorR 3
import_legacy_classifications <- function(path, user_id) {
  old <- utils::read.table(path, header = TRUE, sep = "\t", quote = "", fill = TRUE,
                           colClasses = "character", na.strings = c("", "NA"))
  for (i in seq_len(nrow(old))) {
    add_comment(user_id, variant_key(old$CHROM[i], old$POS[i], old$REF[i], old$ALT[i]),
                classification = old$Classification[i],
                comment = paste0("[PriorR 3 classification by ", old$User[i] %||% "?", "] ", old$Comments[i] %||% ""))
  }
  nrow(old)
}

# ---- "Remember me on this computer" ----
# Only when PriorR runs on the user's own computer (PRIORR_LOCAL, set by launch.R). The
# computer keeps a random token in its local config folder; the shared folder only keeps its
# hash, so copying the shared files is not enough to sign in as somebody else.

remember_enabled <- function() Sys.getenv("PRIORR_LOCAL") == "1"
remember_file <- function() file.path(Sys.getenv("PRIORR_LOGIN_DIR", tools::R_user_dir("priorr", "config")), "login.json")
token_hash <- function(token) sodium::bin2hex(sodium::sha256(charToRaw(token)))

remember_user <- function(user_id) {
  token <- sodium::bin2hex(sodium::random(32))
  modify_user(user_id, function(u) {
    u$remember <- c(u$remember %||% list(), list(list(hash = token_hash(token), created_at = now(),
                                                      computer = Sys.info()[["nodename"]])))
    u
  })
  dir.create(dirname(remember_file()), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(list(username = user_id, token = token), remember_file(), auto_unbox = TRUE)
  Sys.chmod(remember_file(), "0600")
}

# The remembered user of this computer, or NULL
remembered_user <- function() {
  if (!remember_enabled() || !file.exists(remember_file())) return(NULL)
  saved <- tryCatch(jsonlite::fromJSON(remember_file()), error = function(e) NULL)
  if (is.null(saved$username) || is.null(saved$token)) return(NULL)
  a <- read_accounts()
  i <- find_user(a, saved$username)
  if (is.na(i)) return(NULL)
  u <- a$users[[i]]
  hashes <- vapply(u$remember %||% list(), function(r) r$hash %||% "", "")
  if (!isTRUE(u$active) || !token_hash(saved$token) %in% hashes) return(NULL)
  list(id = u$username, username = u$username, full_name = u$full_name, is_admin = isTRUE(u$is_admin))
}

# Sign out: this computer forgets the user
forget_user <- function() {
  if (!file.exists(remember_file())) return(invisible(FALSE))
  saved <- tryCatch(jsonlite::fromJSON(remember_file()), error = function(e) NULL)
  unlink(remember_file())
  if (!is.null(saved$username)) try(modify_user(saved$username, function(u) {
    u$remember <- Filter(function(r) !identical(r$hash, token_hash(saved$token)), u$remember %||% list())
    u
  }), silent = TRUE)
  invisible(TRUE)
}
