# Gene panel library: PanelApp (Genomics England and Australia, downloaded with
# Program_files/scripts/update_panels.R), Orphanet, and the panels saved by the users in the
# shared folder (shared/panels/<name>.txt).
# PanelApp: the lab copy is in Program_files/Panels; a user can download a newer one that only
# they use (shared/panelapp/<user>).

panels_dir <- file.path(data_dir, "Panels")

user_panelapp_dir <- function(username) file.path(shared_dir, "panelapp", username)

# The user of the Shiny session (set at sign in, server.R); NULL outside a session
panel_user <- function() {
  session <- shiny::getDefaultReactiveDomain()
  if (is.null(session)) NULL else session$userData$username
}

# Folder of the PanelApp files used by a user: their own download, or the lab copy
panelapp_dir <- function(username = panel_user()) {
  if (!is.null(username)) {
    dir <- user_panelapp_dir(username)
    if (file.exists(file.path(dir, "panelapp_date.txt"))) return(dir)
  }
  panels_dir
}

panelapp_sources <- c(
  panelapp_genomics_england = "PanelApp England",
  panelapp_australia = "PanelApp Australia")

# Panel files in Program_files/Panels (<file>.tsv) and the name of their source.
# RetiGene: Am J Hum Genet 2025, doi:10.1016/j.ajhg.2025.08.017 (see scripts/import_retigene.R)
panel_file_sources <- c(panelapp_sources, retigene = "RetiGene")

read_panelapp <- function(file, dir = panels_dir) {
  if (!file %in% names(panelapp_sources)) dir <- panels_dir   # RetiGene: always the lab copy
  path <- file.path(dir, paste0(file, ".tsv"))
  if (!file.exists(path)) return(NULL)
  data.table::fread(path, sep = "\t", data.table = FALSE, colClasses = "character")
}

.panel_cache <- new.env()

# All the panels: data.frame(id, source, name, genes = list of green genes, amber = list of amber genes,
# red = list of red genes). PanelApp confidence: 3-4 green, 2 amber, 0-1 red (RetiGene: candidates are amber)
panel_catalog <- function(dir = panelapp_dir()) {
  custom_dir <- file.path(shared_dir, "panels")
  custom_files <- list.files(custom_dir, pattern = "\\.txt$", full.names = TRUE)
  key <- paste(c(dir, file.info(list.files(c(panels_dir, dir), full.names = TRUE))$mtime, file.info(custom_files)$mtime), collapse = "|")
  if (identical(.panel_cache[[dir]]$key, key)) return(.panel_cache[[dir]]$catalog)

  parts <- list()
  for (file in names(panel_file_sources)) {
    pa <- read_panelapp(file, dir)
    if (is.null(pa)) next
    by_panel <- split(pa, pa$panel_id)
    by_panel <- by_panel[order(vapply(by_panel, function(p) p$panel_name[1], ""))]
    parts[[file]] <- data.frame(
      id = paste0(file, ":", names(by_panel)), source = panel_file_sources[[file]],
      name = vapply(by_panel, function(p) sprintf("%s (v%s)", p$panel_name[1], p$version[1]), ""),
      stringsAsFactors = FALSE)
    level <- function(p) suppressWarnings(as.integer(p$confidence))
    parts[[file]]$genes <- lapply(by_panel, function(p) unique(p$gene[level(p) %in% 3:4]))
    parts[[file]]$amber <- lapply(by_panel, function(p) unique(p$gene[level(p) %in% 2]))
    parts[[file]]$red <- lapply(by_panel, function(p) unique(p$gene[level(p) %in% 0:1]))
  }
  orpha <- orphanet_panels("_panel.txt$")
  if (length(orpha)) {
    parts$orphanet <- data.frame(id = paste0("orphanet:", orpha), source = "Orphanet",
                                 name = gsub("_", " ", sub("_panel$", "", orpha)), stringsAsFactors = FALSE)
    parts$orphanet$genes <- lapply(orpha, function(p) read_gene_list(file.path(orphanet_dir, paste0(p, ".txt"))))
    parts$orphanet$amber <- parts$orphanet$red <- lapply(orpha, function(p) character(0))
  }
  # Genes specific of a disease / of a disease class (GLOWgenes, Sci Rep 2026)
  specific <- glow_specific_genes()
  if (nrow(specific)) {
    by_panel <- split(specific$gene, specific$panel)
    parts$glowspecific <- data.frame(id = paste0("glowspecific:", names(by_panel)), source = "GLOWgenes specific genes (Sci Rep 2026)",
                                     name = sub("^(disease|class):", "", names(by_panel)), stringsAsFactors = FALSE)
    parts$glowspecific$name <- ifelse(startsWith(names(by_panel), "class:"), paste(parts$glowspecific$name, "- class-specific genes"),
                                      paste(gsub("_", " ", parts$glowspecific$name), "- disease-specific genes"))
    parts$glowspecific$genes <- unname(lapply(by_panel, unique))
    parts$glowspecific$amber <- parts$glowspecific$red <- lapply(by_panel, function(x) character(0))
  }
  if (length(custom_files)) {
    names <- tools::file_path_sans_ext(basename(custom_files))
    parts$custom <- data.frame(id = paste0("custom:", names), source = "Saved by the users", name = names,
                               stringsAsFactors = FALSE)
    parts$custom$genes <- lapply(custom_files, read_gene_list)
    parts$custom$amber <- parts$custom$red <- lapply(custom_files, function(f) character(0))
  }
  catalog <- if (length(parts)) do.call(rbind, unname(parts)) else
    data.frame(id = character(0), source = character(0), name = character(0))
  rownames(catalog) <- NULL
  .panel_cache[[dir]] <- list(key = key, catalog = catalog)
  catalog
}

# Choices for a selectize grouped by source
# specific: NA = all the panels, FALSE = without, TRUE = only the GLOWgenes specific genes (Benítez et al. 2026)
panel_choices <- function(catalog = panel_catalog(), min_genes = 0, specific = NA) {
  catalog <- catalog[lengths(catalog$genes) >= min_genes, , drop = FALSE]
  is_specific <- startsWith(catalog$id, "glowspecific:")
  if (!is.na(specific)) catalog <- catalog[is_specific == specific, , drop = FALSE]
  if (isTRUE(specific)) catalog$source <- ifelse(startsWith(catalog$id, "glowspecific:class:"),
                                                 "Class-specific genes", "Disease-specific genes")
  if (nrow(catalog) == 0) return(list())
  label <- sprintf("%s [%d genes]", catalog$name, lengths(catalog$genes))
  groups <- unique(catalog$source)
  if (isTRUE(specific)) groups <- c("Disease-specific genes", "Class-specific genes")
  lapply(split(stats::setNames(catalog$id, label), factor(catalog$source, levels = groups)), as.list)
}

# Genes of the selected panels (green genes; also amber ones, and red ones, if asked)
panel_genes <- function(ids, amber = FALSE, catalog = panel_catalog(), red = FALSE) {
  if (!length(ids)) return(NULL)
  sel <- catalog[catalog$id %in% ids, ]
  unique(c(unlist(sel$genes), if (amber || red) unlist(sel$amber), if (red) unlist(sel$red)))
}

# Saves a panel of the user in the shared library
save_custom_panel <- function(name, genes) {
  name <- trimws(name)
  if (!grepl("^[A-Za-z0-9 _.-]+$", name)) stop("Use only letters, numbers, spaces, '-', '_' or '.' in the panel name")
  if (!length(genes)) stop("The panel has no genes")
  dir <- file.path(shared_dir, "panels")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(name, ".txt"))
  if (file.exists(path)) stop("There is already a panel called ", name)
  writeLines(sort(unique(genes)), path)
  invisible(paste0("custom:", name))
}

# Downloads the PanelApp panels (green, amber and red genes) into Program_files/Panels
update_panelapp <- function(progress = function(source, i, n) NULL, dir = panels_dir) {
  apis <- c(panelapp_genomics_england = "https://panelapp.genomicsengland.co.uk/api/v1",
            panelapp_australia = "https://panelapp-aus.org/api/v1")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  failed <- character(0)
  for (file in names(apis)) {
    ids <- character(0)
    url <- paste0(apis[[file]], "/panels/?page=1")
    while (!is.null(url) && !is.na(url)) {
      page <- fetch_json(url)
      if (is.null(page)) stop("PanelApp is not reachable: ", url)
      ids <- c(ids, page$results$id)
      url <- page$`next`
    }
    rows <- vector("list", length(ids))
    for (i in seq_along(ids)) {
      progress(panelapp_sources[[file]], i, length(ids))
      Sys.sleep(0.7)   # PanelApp limits the number of requests per minute
      p <- fetch_json(sprintf("%s/panels/%s/", apis[[file]], ids[i]))
      if (is.null(p)) { failed <- c(failed, paste(panelapp_sources[[file]], ids[i])); next }
      genes <- p$genes
      if (!length(genes) || !NROW(genes)) next
      keep <- !is.na(genes$gene_data$gene_symbol)
      rows[[i]] <- data.frame(panel_id = p$id, panel_name = p$name, version = p$version,
                              gene = genes$gene_data$gene_symbol[keep], confidence = genes$confidence_level[keep],
                              mode_of_inheritance = (genes$mode_of_inheritance %||% rep(NA, nrow(genes)))[keep],
                              stringsAsFactors = FALSE)
    }
    out <- do.call(rbind, rows)
    tmp <- file.path(dir, paste0(file, ".tsv.tmp"))
    data.table::fwrite(out, tmp, sep = "\t")
    file.rename(tmp, file.path(dir, paste0(file, ".tsv")))
  }
  writeLines(format(Sys.Date()), file.path(dir, "panelapp_date.txt"))
  if (length(failed)) warning("Panels that could not be downloaded: ", paste(failed, collapse = ", "))
  invisible(failed)
}

# GET with retries, waiting longer each time (PanelApp answers 429 when there are too many requests)
fetch_json <- function(url, tries = 7) {
  for (i in seq_len(tries)) {
    result <- tryCatch(suppressWarnings(jsonlite::fromJSON(url)), error = function(e) NULL)
    if (!is.null(result)) return(result)
    Sys.sleep(min(2^i, 60))
  }
  NULL
}

# Date of the PanelApp copy used by a user ("lab" or "own")
panelapp_info <- function(username = panel_user()) {
  dir <- panelapp_dir(username)
  date_file <- file.path(dir, "panelapp_date.txt")
  list(date = if (file.exists(date_file)) readLines(date_file, n = 1) else "unknown",
       own = !identical(dir, panels_dir))
}

# Does the PanelApp copy in use have red genes? (copies downloaded before PriorR 4.0 only kept green and amber)
panelapp_has_red <- function(catalog = panel_catalog()) {
  any(lengths(catalog$red[startsWith(catalog$id, "panelapp_")]) > 0)
}
