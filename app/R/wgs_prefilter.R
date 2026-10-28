# WGS prefilter. Large WGS variant files don't fit comfortably in memory, so before loading
# them PriorR:
#   1. Converts the TSV (or TSV.GZ) into a Parquet cache, streaming it in chunks with DuckDB
#      (bounded memory). Done only once per file.
#   2. Applies a prioritization profile over the Parquet cache and writes the candidate
#      variants to a small TSV that is loaded as a normal SNV file.
# It runs in a background R process (Program_files/scripts/wgs_prefilter.R) and reports its
# progress in a JSON file.

sql_str <- function(x) paste0("'", gsub("'", "''", as.character(x), fixed = TRUE), "'")
sql_col <- function(x) paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
sql_num <- function(col) sprintf("TRY_CAST(%s AS DOUBLE)", sql_col(col))

# Cache file of an input file (changes when the file changes)
wgs_cache_path <- function(path, cache_dir) {
  info <- file.info(path)
  key <- paste(normalizePath(path), info$size, as.integer(info$mtime), sep = "|")
  hash <- substr(text_hash(key), 1, 16)
  file.path(cache_dir, paste0(hash, ".parquet"))
}

text_hash <- function(text) {
  f <- tempfile(); on.exit(unlink(f))
  writeLines(text, f, useBytes = TRUE)
  unname(tools::md5sum(f))
}

# SQL of the prioritization profile. A row is kept when it passes the gene filter and either:
#   - the main filter (regions, frequencies, CADD, ClinVar benign), or
#   - a rescue rule (ClinVar pathogenic, or SpliceAI above the threshold with rare frequency).
# Returns list(where, priority).
wgs_profile_sql <- function(cfg, columns) {
  has <- function(col) col %in% columns
  glob <- character(0); main <- character(0); freq <- character(0); rescues <- character(0)

  if (isTRUE(cfg$canonical_only) && has("CANONICAL")) glob <- c(glob, "CANONICAL = 'YES'")
  genes <- unlist(cfg$genes)
  if (length(genes) && has("SYMBOL")) glob <- c(glob, sprintf("SYMBOL IN (%s)", paste(sql_str(genes), collapse = ", ")))
  excluded <- unlist(cfg$exclude_genes)
  if (length(excluded) && has("SYMBOL"))
    glob <- c(glob, sprintf("(SYMBOL IS NULL OR SYMBOL NOT IN (%s))", paste(sql_str(excluded), collapse = ", ")))

  regions <- unlist(cfg$regions)
  if (length(regions) && has("Genomic_region"))
    main <- c(main, sprintf("Genomic_region IN (%s)", paste(sql_str(regions), collapse = ", ")))
  for (col in names(cfg$max_af)) {
    v <- cfg$max_af[[col]]
    if (!is.null(v) && has(col)) freq <- c(freq, sprintf("(%s <= %s OR %s IS NULL)", sql_num(col), as.numeric(v), sql_num(col)))
  }
  if (!is.null(cfg$min_cadd) && has("CADD_PHRED"))
    main <- c(main, sprintf("(%s >= %s OR %s IS NULL)", sql_num("CADD_PHRED"), as.numeric(cfg$min_cadd), sql_num("CADD_PHRED")))
  not_benign <- "coalesce(CLNSIG NOT ILIKE '%benign%', true)"
  if (isTRUE(cfg$exclude_clinvar_benign) && has("CLNSIG")) main <- c(main, not_benign)
  main <- c(main, freq)

  if (isTRUE(cfg$rescue_clinvar_pathogenic) && has("CLNSIG"))
    rescues <- c(rescues, sprintf("(CLNSIG ILIKE '%%pathogenic%%' AND CLNSIG NOT ILIKE '%%conflicting%%' AND %s)", not_benign))
  if (!is.null(cfg$rescue_spliceai_min) && has("SpliceAI_DS_Max"))
    rescues <- c(rescues, sprintf("(%s)", paste(c(sprintf("%s >= %s", sql_num("SpliceAI_DS_Max"), as.numeric(cfg$rescue_spliceai_min)), freq), collapse = " AND ")))

  branches <- c(sprintf("(%s)", if (length(main)) paste(main, collapse = " AND ") else "true"), rescues)
  where <- paste(c(glob, sprintf("(%s)", paste(branches, collapse = " OR "))), collapse = " AND ")

  # ClinVar pathogenic first, then by CADD
  priority <- character(0)
  if (has("CLNSIG")) priority <- c(priority, "(CLNSIG ILIKE '%pathogenic%' AND CLNSIG NOT ILIKE '%conflicting%') DESC")
  if (has("CADD_PHRED")) priority <- c(priority, sprintf("%s DESC NULLS LAST", sql_num("CADD_PHRED")))
  list(where = where, priority = priority)
}

not_a_variant_table <- paste("This file can't be read as a variant table: PriorR opens tab-separated text files",
                             "(.tsv, .txt, also compressed .gz) annotated by the pipeline, with a header line (CHROM, POS, REF, ALT, ...)")

run_wgs_prefilter <- function(input, cfg, output, cache_dir, progress_file, memory_limit = "2GB", threads = 4) {
  start <- Sys.time()
  write_progress <- function(stage, pct = NULL, message = "", ...) {
    tmp <- paste0(progress_file, ".tmp")
    jsonlite::write_json(c(list(stage = stage, pct = pct, message = message), list(...)), tmp,
                         auto_unbox = TRUE, null = "null", digits = NA)
    file.rename(tmp, progress_file)
  }

  tryCatch({
    if (!file.exists(input)) stop("The file doesn't exist: ", input)
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(con, sprintf("SET memory_limit = %s", sql_str(memory_limit)))
    DBI::dbExecute(con, sprintf("SET threads = %d", as.integer(threads)))
    DBI::dbExecute(con, sprintf("SET temp_directory = %s", sql_str(file.path(cache_dir, "duckdb_tmp"))))

    parquet <- wgs_cache_path(input, cache_dir)
    if (!file.exists(parquet)) {
      message <- "Splitting and indexing the WGS file (only the first time)..."
      write_progress("convert", NULL, message)
      options(duckdb.progress_display = function(x) write_progress("convert", x, message))
      DBI::dbExecute(con, "SET enable_progress_bar = true")
      tmp <- paste0(parquet, ".tmp")
      unlink(tmp)
      tryCatch(DBI::dbExecute(con, sprintf(
        "COPY (SELECT * FROM read_csv(%s, delim = '\t', header = true, all_varchar = true, quote = '', escape = ''))
         TO %s (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000)", sql_str(input), sql_str(tmp))),
        error = function(e) {
          unlink(tmp)
          stop(not_a_variant_table, " (", basename(input), ").", call. = FALSE)
        })
      options(duckdb.progress_display = FALSE)
      file.rename(tmp, parquet)
    }

    write_progress("filter", NULL, "Applying the prioritization profile...")
    source_sql <- sprintf("read_parquet(%s)", sql_str(parquet))
    columns <- DBI::dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM %s", source_sql))$column_name
    if (!all(c("CHROM", "POS", "REF", "ALT") %in% columns)) {
      unlink(parquet)
      stop(not_a_variant_table, ": it has no CHROM, POS, REF and ALT columns (", basename(input), ").", call. = FALSE)
    }
    total_rows <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM %s", source_sql))$n
    q <- wgs_profile_sql(cfg, columns)
    DBI::dbExecute(con, sprintf("CREATE TEMP TABLE candidates AS SELECT * FROM %s WHERE %s", source_sql, q$where))
    kept_rows <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM candidates")$n

    # Safety limit (memory): above it nothing is loaded, so a loaded table always has every variant
    # that passes the load filters
    max_rows <- cfg$max_rows
    if (!is.null(max_rows) && kept_rows > max_rows)
      stop(sprintf("%s rows pass the load filters, more than the safety limit of %s rows. Nothing was loaded: use stricter load filters (e.g. the WGS profile, genes or panels), or raise max_rows in priorr.cfg if the computer has enough memory.",
                   format(kept_rows, big.mark = ","), format(as.numeric(max_rows), big.mark = ",", scientific = FALSE)))
    sql <- "SELECT * FROM candidates"

    write_progress("write", NULL, "Writing the candidate variants...")
    df <- DBI::dbGetQuery(con, sql)
    data.table::fwrite(df, paste0(output, ".tmp"), sep = "\t", na = "", quote = FALSE)
    file.rename(paste0(output, ".tmp"), output)

    variants <- if (all(c("CHROM", "POS", "REF", "ALT") %in% names(df))) nrow(unique(df[, c("CHROM", "POS", "REF", "ALT")])) else NA
    write_progress("done", 100, "Done", total_rows = total_rows, kept_rows = kept_rows, loaded_rows = nrow(df),
                   variants = variants,
                   seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1))
    invisible(TRUE)
  }, error = function(e) {
    write_progress("error", NULL, conditionMessage(e))
    invisible(FALSE)
  })
}
