# Projects: groups of variant files (usually one file per sample; a file with several samples
# counts each sample by its <sample>_GT column). While a sample is analysed, each variant says in
# how many other samples of the project it is.
#
# In the shared folder:
#   projects/<id>/project.json        name and files of the project
#   projects/<id>/index/<hash>.parquet variants of a file: CHROM, POS, REF, ALT, sample, gt (het / hom)
#   projects/<id>/index/<hash>.json    the file it comes from (path, size, date) and its samples
#   projects/<id>/indexing.json        progress of the background indexing (scripts/index_project.R)

projects_dir <- function() store_file("projects")
project_dir <- function(id) file.path(projects_dir(), id)

write_json_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE)
  file.rename(tmp, path)
}

read_project <- function(id) {
  path <- file.path(project_dir(id), "project.json")
  if (!file.exists(path)) return(NULL)
  p <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  p$id <- id
  p$files <- p$files %||% list()
  p
}

list_projects <- function() {
  ids <- basename(list.dirs(projects_dir(), recursive = FALSE))
  Filter(Negate(is.null), lapply(stats::setNames(ids, ids), read_project))
}

project_choices <- function() {
  p <- list_projects()
  if (!length(p)) return(character(0))
  stats::setNames(names(p), vapply(p, function(x) x$name, ""))
}

create_project <- function(name, user_id) {
  name <- trimws(name)
  if (!nzchar(name)) stop("Write a name for the project")
  if (tolower(name) %in% tolower(names(project_choices()))) stop("There is already a project called ", name)
  id <- paste0(substr(gsub("[^a-z0-9]+", "_", tolower(name)), 1, 30), "_", substr(text_hash(paste(name, now())), 1, 6))
  write_json_atomic(list(name = name, created_by = user_id, created_at = now(), files = list()),
                    file.path(project_dir(id), "project.json"))
  id
}

delete_project <- function(id) unlink(project_dir(id), recursive = TRUE)

normal_path <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

add_project_files <- function(id, paths, user_id) {
  p <- read_project(id)
  if (is.null(p)) stop("The project doesn't exist")
  have <- vapply(p$files, function(f) normal_path(f$path), "")
  new <- setdiff(unique(normal_path(paths)), have)
  p$files <- c(p$files, lapply(new, function(f) list(path = f, added_by = user_id, added_at = now())))
  write_json_atomic(p[c("name", "created_by", "created_at", "files")], file.path(project_dir(id), "project.json"))
  length(new)
}

remove_project_file <- function(id, path) {
  p <- read_project(id)
  p$files <- Filter(function(f) normal_path(f$path) != normal_path(path), p$files)
  write_json_atomic(p[c("name", "created_by", "created_at", "files")], file.path(project_dir(id), "project.json"))
  unlink(paste0(index_base(id, path), c(".parquet", ".json")))
}

# Index of one file of the project (without extension)
index_base <- function(id, path) file.path(project_dir(id), "index", substr(text_hash(normal_path(path)), 1, 16))

# Is the index of a file up to date (same size and date of the file)? NULL when not indexed.
index_info <- function(id, path) {
  meta <- paste0(index_base(id, path), ".json")
  if (!file.exists(meta) || !file.exists(paste0(index_base(id, path), ".parquet"))) return(NULL)
  m <- tryCatch(jsonlite::fromJSON(meta, simplifyVector = TRUE), error = function(e) NULL)
  info <- file.info(path)
  if (is.null(m) || is.na(info$size) || !identical(as.numeric(m$size), as.numeric(info$size)) ||
      !identical(as.numeric(m$mtime), as.numeric(as.integer(info$mtime)))) return(NULL)
  m
}

# Files of a project with their state: indexed (and samples) or pending
project_status <- function(id) {
  p <- read_project(id)
  if (is.null(p) || !length(p$files)) return(data.frame(path = character(0), indexed = logical(0), samples = character(0)))
  paths <- vapply(p$files, `[[`, "", "path")
  infos <- lapply(paths, function(f) index_info(id, f))
  data.frame(path = paths, exists = file.exists(paths), indexed = !vapply(infos, is.null, TRUE),
             samples = vapply(infos, function(m) paste(unlist(m$samples), collapse = ", "), ""),
             stringsAsFactors = FALSE)
}

# Samples of the indexed files of a project
project_samples <- function(id) {
  st <- project_status(id)
  unique(unlist(strsplit(st$samples[st$indexed], ", ", fixed = TRUE)))
}

# SQL of the genotype of a sample column: het (one alternative allele), hom (all alternative,
# also hemizygous "1"), NULL (no alternative allele)
gt_sql <- function(col) {
  g <- sprintf("TRIM(CAST(%s AS VARCHAR))", sql_col(col))
  sprintf("CASE WHEN NOT regexp_matches(%s, '[1-9]') THEN NULL
                WHEN regexp_matches(%s, '(^|[/|])[0.]([/|]|$)') THEN 'het' ELSE 'hom' END", g, g)
}

# Index of a variant file: the variants (CHROM, POS, REF, ALT) of each sample with its genotype.
# The cache of PriorR (Parquet of a big file opened before) is used when it is there. Without
# <sample>_GT columns the file is one sample named as the file.
index_variant_file <- function(path, base, cache_dir = wgs_cache_dir, memory_limit = "2GB", threads = 4) {
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, sprintf("SET memory_limit = %s", sql_str(memory_limit)))
  DBI::dbExecute(con, sprintf("SET threads = %d", as.integer(threads)))
  cached <- wgs_cache_path(path, cache_dir)
  source_sql <- if (file.exists(cached)) sprintf("read_parquet(%s)", sql_str(cached))
                else sprintf("read_csv(%s, delim = '\t', header = true, all_varchar = true, quote = '', escape = '')", sql_str(path))
  columns <- DBI::dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM %s", source_sql))$column_name
  if (!all(c("CHROM", "POS", "REF", "ALT") %in% columns))
    stop("It has no CHROM, POS, REF and ALT columns: ", basename(path), call. = FALSE)
  samples <- sample_names(columns)
  # POS as an integer in text (also when a file writes 1e+08), like key_text
  key <- paste("CAST(CHROM AS VARCHAR) AS CHROM,",
               "COALESCE(CAST(CAST(TRY_CAST(POS AS DOUBLE) AS BIGINT) AS VARCHAR), CAST(POS AS VARCHAR)) AS POS,",
               "CAST(REF AS VARCHAR) AS REF, CAST(ALT AS VARCHAR) AS ALT")
  parts <- if (length(samples)) vapply(samples, function(s) sprintf(
      "SELECT DISTINCT %s, %s AS sample, %s AS gt FROM %s WHERE (%s) IS NOT NULL",
      key, sql_str(s), gt_sql(paste0(s, "_GT")), source_sql, gt_sql(paste0(s, "_GT"))), "")
    else sprintf("SELECT DISTINCT %s, %s AS sample, NULL AS gt FROM %s", key,
                 sql_str(sub("\\.(tsv|txt|csv)(\\.gz)?$", "", basename(path))), source_sql)
  if (!length(samples)) samples <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", basename(path))
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(base, ".parquet.tmp", Sys.getpid())
  DBI::dbExecute(con, sprintf("COPY (%s ORDER BY CHROM, POS) TO %s (FORMAT parquet, COMPRESSION zstd)",
                              paste(parts, collapse = " UNION ALL "), sql_str(tmp)))
  file.rename(tmp, paste0(base, ".parquet"))
  info <- file.info(path)
  write_json_atomic(list(path = normal_path(path), size = info$size, mtime = as.integer(info$mtime),
                         samples = as.list(samples), indexed_at = now()), paste0(base, ".json"))
  invisible(samples)
}

# Indexes the files of a project that are new or changed (scripts/index_project.R runs it in the
# background); progress in indexing.json
index_project <- function(id, cache_dir = wgs_cache_dir, memory_limit = "2GB", threads = 4) {
  progress <- file.path(project_dir(id), "indexing.json")
  st <- project_status(id)
  todo <- st$path[!st$indexed & st$exists]
  errors <- character(0)
  for (i in seq_along(todo)) {
    write_json_atomic(list(stage = "running", done = i - 1, total = length(todo), file = basename(todo[i]),
                           pid = Sys.getpid(), time = now()), progress)
    result <- tryCatch(index_variant_file(todo[i], index_base(id, todo[i]), cache_dir, memory_limit, threads),
                       error = function(e) e)
    if (inherits(result, "error")) errors <- c(errors, paste0(basename(todo[i]), ": ", conditionMessage(result)))
  }
  write_json_atomic(list(stage = "done", done = length(todo), total = length(todo), errors = as.list(errors), time = now()), progress)
  invisible(length(errors) == 0)
}

# Background indexing in progress? (the process is still alive)
project_indexing <- function(id) {
  path <- file.path(project_dir(id), "indexing.json")
  if (!file.exists(path)) return(NULL)
  s <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  if (identical(s$stage, "running")) {
    if (!process_alive(s$pid)) s$stage <- "stopped"
  }
  s
}

# (signal 0 only checks on Unix; on Windows pskill would stop the process, so tasklist)
process_alive <- function(pid) {
  if (is.null(pid) || is.na(pid)) return(FALSE)
  if (.Platform$OS.type == "windows") {
    out <- tryCatch(system2("tasklist", c("/FI", shQuote(sprintf("PID eq %d", as.integer(pid))), "/NH"), stdout = TRUE),
                    error = function(e) character(0))
    return(any(grepl(sprintf("\\b%d\\b", as.integer(pid)), out)))
  }
  isTRUE(tools::pskill(as.integer(pid), 0L))
}

start_project_indexing <- function(id) {
  s <- project_indexing(id)
  if (!is.null(s) && identical(s$stage, "running")) return(FALSE)
  script <- file.path(data_dir, "scripts", "index_project.R")
  log <- file.path(project_dir(id), "indexing.log")
  system2(file.path(R.home("bin"), "Rscript"),
          c(shQuote(normalizePath(script)), shQuote(id)), stdout = log, stderr = log, wait = FALSE)
  TRUE
}

# Text form of the key columns, as in the index (POS without scientific notation)
key_text <- function(x) if (is.numeric(x)) sprintf("%.0f", x) else as.character(x)

# For the variants of `df`: in how many samples of the project (not in `exclude`) each one is,
# and in how many of them homozygous. Returns list(counts = data.frame(key, n, hom), total = samples).
project_counts <- function(id, df, exclude = character(0)) {
  st <- project_status(id)
  files <- vapply(st$path[st$indexed], function(f) paste0(index_base(id, f), ".parquet"), "")
  samples <- setdiff(project_samples(id), exclude)
  empty <- list(counts = data.frame(key = character(0), n = integer(0), hom = integer(0)), total = length(samples))
  if (!length(files) || !nrow(df)) return(empty)
  keys <- unique(data.frame(CHROM = key_text(df$CHROM), POS = key_text(df$POS), REF = key_text(df$REF),
                            ALT = key_text(df$ALT), stringsAsFactors = FALSE))
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  duckdb::duckdb_register(con, "k", keys)
  out <- DBI::dbGetQuery(con, sprintf(
    "SELECT i.CHROM, i.POS, i.REF, i.ALT, count(DISTINCT i.sample) AS n,
            count(DISTINCT i.sample) FILTER (WHERE i.gt = 'hom') AS hom
     FROM read_parquet([%s]) i JOIN k USING (CHROM, POS, REF, ALT)
     %s GROUP BY ALL",
    paste(sql_str(files), collapse = ", "),
    if (length(exclude)) sprintf("WHERE i.sample NOT IN (%s)", paste(sql_str(exclude), collapse = ", ")) else ""))
  list(counts = data.frame(key = variant_key(out$CHROM, out$POS, out$REF, out$ALT), n = as.integer(out$n),
                           hom = as.integer(out$hom), stringsAsFactors = FALSE),
       total = length(samples))
}

# Samples of the project with one variant, their genotype and file
project_variant_samples <- function(id, chrom, pos, ref, alt) {
  st <- project_status(id)
  st <- st[st$indexed, , drop = FALSE]
  if (!nrow(st)) return(data.frame(sample = character(0), gt = character(0), file = character(0)))
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  sql <- paste(vapply(st$path, function(f) sprintf(
    "SELECT sample, gt, %s AS file FROM read_parquet(%s) WHERE CHROM = %s AND POS = %s AND REF = %s AND ALT = %s",
    sql_str(f), sql_str(paste0(index_base(id, f), ".parquet")), sql_str(key_text(chrom)), sql_str(key_text(pos)),
    sql_str(key_text(ref)), sql_str(key_text(alt))), ""), collapse = " UNION ALL ")
  DBI::dbGetQuery(con, paste("SELECT DISTINCT * FROM (", sql, ") ORDER BY sample"))
}
