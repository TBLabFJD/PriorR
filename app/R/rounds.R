# Analysis rounds of a file: the variants shown with some filters are saved as a round, so the
# next filters can hide them (only the new variants are seen). Per user and file, in the shared
# folder (rounds/<user>/<hash of the file path>.json).

rounds_path <- function(user_id, file) {
  dir <- store_file("rounds", sub("\\.jsonl$", "", user_file_name(user_id)))
  file.path(dir, paste0(substr(text_hash(normalizePath(file, mustWork = FALSE)), 1, 16), ".json"))
}

list_rounds <- function(user_id, file) {
  path <- rounds_path(user_id, file)
  if (!file.exists(path)) return(list())
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE)$rounds, error = function(e) list())
}

write_rounds <- function(user_id, file, rounds) {
  path <- rounds_path(user_id, file)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(list(file = file, rounds = rounds), tmp, auto_unbox = TRUE, null = "null")
  file.rename(tmp, path)
}

add_round <- function(user_id, file, keys, name = "", summary = "") {
  rounds <- list_rounds(user_id, file)
  keys <- unique(keys)
  if (!length(keys)) stop("There are no variants to save")
  if (!nzchar(trimws(name))) name <- paste("Round", length(rounds) + 1)
  rounds[[length(rounds) + 1]] <- list(name = trimws(name), created_at = now(), summary = summary, keys = as.list(keys))
  write_rounds(user_id, file, rounds)
  invisible(length(rounds))
}

delete_round <- function(user_id, file, i) {
  rounds <- list_rounds(user_id, file)
  if (i >= 1 && i <= length(rounds)) write_rounds(user_id, file, rounds[-i])
}

# Variant keys of all the rounds of a file
round_keys <- function(user_id, file) unique(unlist(lapply(list_rounds(user_id, file), `[[`, "keys")))

# Name of a round from the sample and the filters: "P001 · Retinal disorders (v4.2), RP +1 · AF <= 0.01 · CADD >= 20"
round_name <- function(sample = NA, panels = character(0), hpo = 0, genes = character(0), max_af = NULL, min_cadd = 0) {
  short <- function(x, n = 2) if (length(x) > n) paste0(paste(x[1:n], collapse = ", "), " +", length(x) - n) else paste(x, collapse = ", ")
  parts <- c(if (!is.na(sample) && nzchar(sample %||% "")) sample,
             if (length(panels)) short(panels),
             if (hpo > 0) sprintf("%d HPO", hpo),
             if (length(genes)) short(genes, 3),
             if (!is.null(max_af)) sprintf("AF \u2264 %s", format(max_af)),
             if (!is.na(min_cadd) && min_cadd > 0) sprintf("CADD \u2265 %s", format(min_cadd)))
  if (length(parts) <= 1) parts <- c(parts, "no filters")
  paste(parts, collapse = " \u00b7 ")
}
