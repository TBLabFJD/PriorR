# GLOWgenes ranking view: rare variants with evidence of pathogenicity, ordered by the GLOWgenes
# rank of their gene (0 = genes of the disease panel, then the closest genes). The rank comes from
# the GLOWgenes column of the file or from the ranking of an Orphanet panel. The filters are in
# Program_files/glowgenes.json.

glow_settings <- function(path = file.path(data_dir, "glowgenes.json")) jsonlite::fromJSON(path, simplifyVector = FALSE)

has_glow_column <- function(df) "GLOWgenes" %in% names(df) && any(!is.na(suppressWarnings(as.numeric(df$GLOWgenes))))

# ---- Rankings available ----
# "orphanet:<panel>"  GLOWgenes rankings of Orphanet panels (Paneles_Orphanet/*_panel_glow*.txt)
# "panelapp:<panel>"  GLOWgenes rankings of 209 PanelApp panels (Sci Rep 2026, Program_files/GLOWgenes)
# "run:<name>"        rankings computed by the lab with GLOWgenes (shared/glowgenes/<name>.tsv)

glow_dir <- file.path(data_dir, "GLOWgenes")
glow_panelapp_file <- file.path(glow_dir, "glowgenes_panelapp_ranks.tsv.gz")
glow_runs_dir <- function() file.path(shared_dir, "glowgenes")

# Runs of GLOWgenes in progress or failed (glowgenes/work/<name>/status.json)
glow_run_statuses <- function() {
  files <- Sys.glob(file.path(glow_runs_dir(), "work", "*", "status.json"))
  rows <- lapply(files, function(f) {
    st <- tryCatch(jsonlite::fromJSON(f), error = function(e) list())
    data.frame(name = basename(dirname(f)), stage = st$stage %||% "running", message = st$message %||% "",
               user = st$user %||% "", time = substr(st$time %||% "", 1, 16), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, c(list(data.frame(name = character(0), stage = character(0), message = character(0),
                                          user = character(0), time = character(0))), rows))
  out[out$stage != "done", , drop = FALSE]
}

glow_panelapp_panels <- function() {
  if (!file.exists(glow_panelapp_file)) return(data.frame(panel = character(0), class = character(0)))
  header <- strsplit(readLines(gzcon(file(glow_panelapp_file, "rb")), n = 1), "\t")[[1]][-1]
  classes <- read_tsv_gz(file.path(glow_dir, "glowgenes_panelapp_panels.tsv.gz"))
  data.frame(panel = header, class = classes$Disease_Class[match(header, classes$PanelApp_Panel)] %||% NA, stringsAsFactors = FALSE)
}

glow_rank_choices <- function() {
  label <- function(x) gsub("_", " ", x)
  orpha <- orphanet_panels()
  pa <- glow_panelapp_panels()
  runs <- tools::file_path_sans_ext(list.files(glow_runs_dir(), pattern = "\\.tsv$"))
  c(list(`Orphanet` = as.list(stats::setNames(paste0("orphanet:", orpha), label(sub("_panel_glow$", "", orpha))))),
    if (nrow(pa)) list(`PanelApp (GLOWgenes, Sci Rep 2026)` = as.list(stats::setNames(paste0("panelapp:", pa$panel),
      ifelse(is.na(pa$class), label(pa$panel), sprintf("%s (%s)", label(pa$panel), pa$class))))),
    if (length(runs)) list(`Computed by the lab` = as.list(stats::setNames(paste0("run:", runs), runs))))
}

.glow_cache <- new.env()

# Gene ranks of one ranking: data.frame(SYMBOL, rank)
glow_ranking <- function(id) {
  if (!is.null(.glow_cache[[id]])) return(.glow_cache[[id]])
  source <- sub(":.*$", "", id)
  name <- sub("^[^:]*:", "", id)
  r <- switch(source,
    orphanet = { x <- orphanet_rank(name); data.frame(SYMBOL = x$SYMBOL, rank = x$Orphanet_GLOW_rank) },
    panelapp = {
      con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE)); on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
      x <- DBI::dbGetQuery(con, sprintf("SELECT SYMBOL, %s AS rank FROM read_csv(%s, delim = '\t', header = true, all_varchar = true)",
                                        sql_col(name), sql_str(glow_panelapp_file)))
      data.frame(SYMBOL = x$SYMBOL, rank = suppressWarnings(as.numeric(x$rank)))
    },
    run = data.table::fread(file.path(glow_runs_dir(), paste0(name, ".tsv")), data.table = FALSE)[, c("SYMBOL", "rank")],
    stop("Unknown ranking ", id))
  .glow_cache[[id]] <- r
  r
}

# Rank of each row: from the selected rankings (the best one when several), or the GLOWgenes column
glow_rank <- function(DF, sources = NULL) {
  if (length(sources)) {
    ranks <- do.call(rbind, lapply(sources, glow_ranking))
    best <- stats::aggregate(rank ~ SYMBOL, data = ranks[!is.na(ranks$rank), ], FUN = min)
    return(best$rank[match(DF$SYMBOL, best$SYMBOL)])
  }
  if ("GLOWgenes" %in% names(DF)) suppressWarnings(as.numeric(DF$GLOWgenes)) else rep(NA_real_, nrow(DF))
}

# ---- Running GLOWgenes (Docker image tblabfjd/glowgenes) ----

glow_docker_image <- "tblabfjd/glowgenes"

# Smallest gene list for a GLOWgenes run
glow_min_genes <- 10

docker_available <- function() !nzchar(docker_problem())

# Why Docker can't be used ("" when it can)
docker_problem <- function() {
  docker <- Sys.which("docker")
  if (!nzchar(docker)) return("Docker is not installed on this computer: install Docker Desktop (Windows, macOS) or Docker Engine (Linux), then try again.")
  out <- suppressWarnings(system2(docker, "info", stdout = TRUE, stderr = TRUE))
  if (is.null(attr(out, "status"))) return("")
  if (any(grepl("permission denied", out, ignore.case = TRUE)))
    return("Your user can't use Docker: on Linux run 'sudo usermod -aG docker $USER', log out and log in again.")
  "Docker is installed but not running: open Docker Desktop (or start the Docker service) and try again."
}

# Ranking of a GLOWgenes output: panel genes rank 0, then GLOWgenes_prioritization_Random.txt (1...)
glow_output_ranking <- function(prioritization_file, panel_genes) {
  p <- data.table::fread(prioritization_file, header = FALSE, data.table = FALSE)
  ranked <- data.frame(SYMBOL = as.character(p[[1]]), rank = as.numeric(p[[ncol(p)]]), stringsAsFactors = FALSE)
  rbind(data.frame(SYMBOL = unique(panel_genes), rank = 0, stringsAsFactors = FALSE),
        ranked[!ranked$SYMBOL %in% panel_genes, ])
}

glow_view <- function(DF, sources = NULL, s = glow_settings()) {
  num <- function(col) if (col %in% names(DF)) suppressWarnings(as.numeric(DF[[col]])) else rep(NA_real_, nrow(DF))
  keep <- rep(TRUE, nrow(DF))
  for (col in names(s$max_af)) if (!is.null(s$max_af[[col]])) keep <- keep & (num(col) <= s$max_af[[col]] | is.na(num(col)))
  clnsig <- if ("CLNSIG" %in% names(DF)) as.character(DF$CLNSIG) else rep(NA_character_, nrow(DF))
  if (isTRUE(s$exclude_clinvar_benign)) keep <- keep & !grepl("benign", clnsig, ignore.case = TRUE)

  e <- s$evidence
  priority <- priority_class(DF)
  evidence <- rep(FALSE, nrow(DF))
  if (isTRUE(e$clinvar_pathogenic)) evidence <- evidence | (grepl("pathogenic", clnsig, ignore.case = TRUE) & !grepl("conflicting", clnsig, ignore.case = TRUE))
  if (isTRUE(e$high_impact) && "Consequence" %in% names(DF))
    evidence <- evidence | (if ("IMPACT" %in% names(DF)) DF$IMPACT %in% "HIGH" | consequence_impact(DF$Consequence) == "HIGH"
                            else consequence_impact(DF$Consequence) == "HIGH")
  if (!is.null(e$min_cadd)) evidence <- evidence | num("CADD_PHRED") >= e$min_cadd
  if (!is.null(e$min_revel)) evidence <- evidence | num("REVELScore") >= e$min_revel
  if (!is.null(e$min_spliceai)) evidence <- evidence | num("SpliceAI_DS_Max") >= e$min_spliceai
  if (isTRUE(e$priorr_priority)) evidence <- evidence | !is.na(priority)

  DF$GLOW_rank <- glow_rank(DF, sources)
  keep <- keep & (evidence %in% TRUE)
  if (!is.null(s$max_rank)) keep <- keep & DF$GLOW_rank <= s$max_rank
  DF <- DF[keep %in% TRUE, , drop = FALSE]
  priority <- priority[keep %in% TRUE]

  # One row per variant: the canonical transcript when there is one
  if (isTRUE(s$one_row_per_variant) && nrow(DF)) {
    canonical_first <- order(!(DF$CANONICAL %in% "YES"))
    DF <- DF[canonical_first, , drop = FALSE]
    priority <- priority[canonical_first]
    first <- !duplicated(row_variant_ids(DF))
    DF <- DF[first, , drop = FALSE]
    priority <- priority[first]
  }

  cadd <- if ("CADD_PHRED" %in% names(DF)) suppressWarnings(as.numeric(DF$CADD_PHRED)) else rep(NA_real_, nrow(DF))
  o <- order(DF$GLOW_rank, match(priority, c("red", "orange"), nomatch = 3L), -cadd, na.last = TRUE)
  DF[o, , drop = FALSE]
}

# ---- Genes specific of diseases / disease classes and gene-disease specificity score (SGDS) ----
# From Benitez et al., Sci Rep 2026 (Program_files/GLOWgenes, see its README)

glow_specific_genes <- function() {
  read <- function(f) { path <- file.path(glow_dir, f); if (file.exists(path)) read_tsv_gz(path) else NULL }
  d <- read("specific_genes_disease.tsv.gz")
  c <- read("specific_genes_class.tsv.gz")
  rbind(if (!is.null(d)) data.frame(gene = trimws(d$Gene), panel = paste0("disease:", trimws(d$Gene_Panel))),
        if (!is.null(c)) data.frame(gene = trimws(c$Gene), panel = paste0("class:", trimws(c$Disease_Class))),
        data.frame(gene = character(0), panel = character(0)))
}

.sgds <- new.env()
sgds_table <- function() {
  if (is.null(.sgds$table)) {
    path <- file.path(glow_dir, "sgds.tsv.gz")
    .sgds$table <- if (file.exists(path)) read_tsv_gz(path) else data.frame(Symbols = character(0), SGDS = numeric(0))
  }
  .sgds$table
}

# Adds the SGDS of the gene of each row (column SGDS)
add_sgds <- function(df) {
  if (!"SYMBOL" %in% names(df) || "SGDS" %in% names(df)) return(df)
  t <- sgds_table()
  df$SGDS <- round(t$SGDS[match(df$SYMBOL, t$Symbols)], 3)
  df
}

# Best GLOWgenes rank of the genes of each CNV (Gene_name "A;B;C")
cnv_glow_rank <- function(gene_field, sources) {
  ranks <- do.call(rbind, lapply(sources, glow_ranking))
  ranks <- ranks[!is.na(ranks$rank), ]
  best <- tapply(ranks$rank, ranks$SYMBOL, min)
  genes <- strsplit(as.character(gene_field), "[,;|&]")
  row <- rep(seq_along(genes), lengths(genes))
  r <- best[match(trimws(unlist(genes)), names(best))]
  out <- rep(NA_real_, length(genes))
  ok <- !is.na(r)
  if (any(ok)) {
    m <- tapply(r[ok], row[ok], min)
    out[as.integer(names(m))] <- m
  }
  out
}
