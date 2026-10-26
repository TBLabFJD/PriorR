# Filtering profiles. A profile has:
#   load:    filters applied while reading the file (see wgs_prefilter.R), to load only candidates
#   filters: values of the SNV sidebar filters (named as the inputs of mod_snv.R)
#   columns: displayed columns (NULL = the default ones)
# Built-in profiles are in Program_files/profiles; each user saves their own ones in
# shared/profiles/<user>.json (only that user writes the file).

profiles_dir <- file.path(data_dir, "profiles")

# Sidebar filters that don't filter anything
default_filters <- function() {
  list(type = NULL, region = NULL, consequence = NULL, canonical = FALSE, exclude_predicted = FALSE,
       max_gnomadg_af = 1, max_gnomadg_popmax = 1, max_gnomade_af = 1, max_fjd_af = 1, min_cadd = 0,
       exclude_benign = FALSE, inheritance = "", panels = list(), amber = FALSE, red = FALSE, excluded_panels = list(), phenotypes = list(),
       orphanet = list(), glow = -1, color_priority = FALSE, prioritize = FALSE)
}

default_load <- function() {
  list(regions = list(), canonical_only = FALSE, max_af = list(), min_cadd = NULL, exclude_clinvar_benign = FALSE,
       rescue_clinvar_pathogenic = FALSE, rescue_spliceai_min = NULL, genes = list(), panels = list(),
       amber = FALSE, red = FALSE, excluded_panels = list())
}

# Values of `x` over the defaults (one level: lists such as the regions are replaced, not merged)
with_defaults <- function(defaults, x) {
  for (n in names(x)) defaults[n] <- list(x[[n]])
  defaults
}

complete_profile <- function(p) {
  p$load <- with_defaults(default_load(), p$load)
  p$filters <- with_defaults(default_filters(), p$filters)
  p
}

builtin_profiles <- function() {
  files <- sort(list.files(profiles_dir, pattern = "\\.json$", full.names = TRUE))
  profiles <- lapply(files, function(f) complete_profile(jsonlite::fromJSON(f, simplifyVector = FALSE)))
  stats::setNames(profiles, vapply(profiles, `[[`, "", "name"))
}

user_profiles_file <- function(user_id) file.path(shared_dir, "profiles", sub("\\.jsonl$", ".json", user_file_name(user_id)))

user_profiles <- function(user_id) {
  f <- user_profiles_file(user_id)
  if (!file.exists(f)) return(list())
  profiles <- lapply(jsonlite::fromJSON(f, simplifyVector = FALSE), complete_profile)
  stats::setNames(profiles, vapply(profiles, `[[`, "", "name"))
}

write_user_profiles <- function(user_id, profiles) {
  f <- user_profiles_file(user_id)
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(f, ".tmp")
  jsonlite::write_json(unname(profiles), tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  file.rename(tmp, f)
}

# Saves (or replaces) a profile of the user
save_profile <- function(user_id, profile) {
  name <- trimws(profile$name %||% "")
  if (!nzchar(name)) stop("Give the profile a name")
  if (name %in% names(builtin_profiles())) stop("'", name, "' is a built-in profile: choose another name")
  profile$name <- name
  profiles <- user_profiles(user_id)
  profiles[[name]] <- profile
  write_user_profiles(user_id, profiles)
  invisible(name)
}

delete_profile <- function(user_id, name) {
  profiles <- user_profiles(user_id)
  profiles[[name]] <- NULL
  write_user_profiles(user_id, profiles)
}

# Built-in and own profiles; ids "builtin:<name>" and "user:<name>"
all_profiles <- function(user_id = NULL) {
  b <- builtin_profiles()
  u <- if (is.null(user_id)) list() else user_profiles(user_id)
  c(stats::setNames(b, sprintf("builtin:%s", names(b))), stats::setNames(u, sprintf("user:%s", names(u))))
}

profile_choices <- function(user_id = NULL) {
  p <- all_profiles(user_id)
  ids <- names(p)
  labels <- vapply(p, `[[`, "", "name")
  mine <- startsWith(ids, "user:")
  c(stats::setNames(ids[!mine], labels[!mine]), stats::setNames(ids[mine], sprintf("%s (mine)", labels[mine])))
}

json_num <- function(x) if (is.null(x)) NA else x
