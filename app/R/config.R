# Configuration. Each setting is taken from (first found):
#   1. an environment variable (PRIORR_SHARED_DIR...), for advanced uses
#   2. priorr.cfg of the user of this computer (written by the first-run setup, see mod_setup.R)
#   3. priorr.cfg next to PriorR (launch.R): the configuration of the whole lab
#   4. the defaults
# priorr.cfg is a text file with "key: value" lines (see priorr.cfg.example).
# Other settings: genome_build (links.R).

priorr_version <- "4.2.1"
priorr_name <- paste0("PriorR v", priorr_version)

priorr_home <- normalizePath("..", mustWork = FALSE)        # the app runs in <PriorR>/app
user_config_file <- function() Sys.getenv("PRIORR_USER_CONFIG", file.path(tools::R_user_dir("priorr", "config"), "priorr.cfg"))
lab_config_file <- function() file.path(priorr_home, "priorr.cfg")

# "key: value" lines; lines starting with # and empty lines are ignored
read_config <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(list())
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#") & grepl(":", lines, fixed = TRUE)]
  keys <- trimws(sub(":.*$", "", lines))
  stats::setNames(as.list(trimws(sub("^[^:]*:", "", lines))), keys)
}

# Folders of the lab configuration can be relative to the PriorR folder ("shared_folder: shared"),
# so a copy of the PriorR folder with its shared folder works anywhere
path_settings <- c("shared_folder", "cache_folder", "alphamissense_file")
is_absolute_path <- function(x) grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", x)

setting <- function(key, env, default) {
  if (nzchar(Sys.getenv(env))) return(Sys.getenv(env))
  configs <- list(user = read_config(user_config_file()), lab = read_config(lab_config_file()))
  for (origin in names(configs)) {
    value <- configs[[origin]][[key]]
    if (is.null(value) || !nzchar(value)) next
    if (origin == "lab" && key %in% path_settings && !is_absolute_path(value)) value <- file.path(priorr_home, value)
    return(path.expand(value))
  }
  default
}

# Data shipped with PriorR (panels, HPO, profiles...)
data_dir <- Sys.getenv("PRIORR_DIR", "Program_files")
orphanet_dir <- file.path(data_dir, "Paneles_Orphanet")
wgs_script <- file.path(data_dir, "scripts", "wgs_prefilter.R")

# Settings that can change (the first-run setup writes the user's priorr.cfg and reloads them)
shared_dir <- wgs_cache_dir <- wgs_memory_limit <- wgs_threads <- data_folders <- NULL
load_settings <- function() {
  # Shared folder of the lab: users, groups, comments, profiles, saved sessions (see store.R)
  shared_dir <<- setting("shared_folder", "PRIORR_SHARED_DIR", file.path(priorr_home, "shared"))
  # Big variant files are cached in the local disk (faster)
  wgs_cache_dir <<- setting("cache_folder", "PRIORR_CACHE_DIR", tools::R_user_dir("priorr", "cache"))
  wgs_memory_limit <<- setting("memory_limit", "PRIORR_WGS_MEMORY", "2GB")
  wgs_threads <<- setting("threads", "PRIORR_WGS_THREADS", "4")
  data_folders <<- setting("data_folders", "PRIORR_DATA_DIRS", "")
  # Safety limit of rows loaded (memory): above it nothing is loaded
  max_rows_setting <<- as.numeric(setting("max_rows", "PRIORR_MAX_ROWS", "2000000"))
  # AlphaMissense predictions (scripts/get_alphamissense.R), only used when the file doesn't have them
  alphamissense_file <<- setting("alphamissense_file", "PRIORR_ALPHAMISSENSE",
                                 file.path(shared_dir, "annotations", "AlphaMissense_hg38.parquet"))
  invisible(TRUE)
}
load_settings()

# Is the shared folder configured somewhere (or PriorR runs as a server, not on a user's computer)?
shared_folder_configured <- function() {
  Sys.getenv("PRIORR_LOCAL") != "1" || nzchar(Sys.getenv("PRIORR_SHARED_DIR")) ||
    !is.null(read_config(user_config_file())$shared_folder) || !is.null(read_config(lab_config_file())$shared_folder)
}

# Writes settings in the user's priorr.cfg (keeping the other ones)
save_user_config <- function(values) {
  current <- read_config(user_config_file())
  current[names(values)] <- values
  dir.create(dirname(user_config_file()), showWarnings = FALSE, recursive = TRUE)
  writeLines(c("# PriorR configuration of this computer (see the Help of PriorR)",
               sprintf("%s: %s", names(current), unlist(current))), user_config_file())
  load_settings()
}

# Folders where "Open file" can browse, as "Name=path;Name2=path2".
# By default: the home folder, the PriorR shared folder and the disks of the computer.
parse_roots <- function(value) {
  pairs <- strsplit(trimws(strsplit(value, ";", fixed = TRUE)[[1]]), "=", fixed = TRUE)
  stats::setNames(vapply(pairs, function(p) path.expand(trimws(p[2])), ""), vapply(pairs, function(p) trimws(p[1]), ""))
}
file_roots <- function() {
  custom <- if (nzchar(data_folders)) parse_roots(data_folders) else character(0)
  c(custom, Home = path.expand("~"), `Shared folder` = shared_dir, shinyFiles::getVolumes()())
}

# Small uploads (panels, BED, PED) up to 2.5 GB; big files are opened from the disk
options(shiny.maxRequestSize = 2500 * 1024^2)

genomic_regions <- c("3UTR", "5UTR", "DOWNSTREAM", "EXONIC", "INTRONIC", "ncRNA", "regulatory", "SPLICING", "UPSTREAM")
inheritance_modes <- c("De novo", "Dominant", "Recessive-Homozygous", "Recessive-Compound heterozygous", "X-linked")
classifications <- c("Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic")

# Columns shown by default when a SNV file is loaded (sample columns are added automatically)
default_snv_columns <- c(
  "CHROM", "POS", "REF", "ALT", "Location", "SYMBOL", "Gene_full_name", "VARIANT_CLASS",
  "Existing_variation", "Genomic_region", "CANONICAL", "Feature", "BIOTYPE", "Consequence",
  "HGVSc", "HGVSp", "CLNSIG", "OMIM_phenotype", "Orphanet_disorder", "CADD_PHRED",
  "gnomADg_AF", "gnomADg_AF_popmax", "gnomADg_popmax", "gnomADe_AF", "FJD_MAF_AF", "CSVS_AF",
  "MutScore", "ada_score", "rf_score", "N_Pathogenic_pred", "N_Benign_pred", "N_predictions",
  "kaviar_AF", "ExACpLI", "gnomAD_LOEUF", "AlphaMissense", "SGDS")

default_cnv_columns <- c(
  "ACMG_class", "Samples_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type",
  "N_PROGRAMS", "Annotation_mode", "Gene_name", "Gene_count", "Tx", "Frameshift", "Exon_count",
  "Location", "Location2", "P_gain_phen", "P_gain_hpo", "P_loss_phen", "P_loss_hpo",
  "P_snvindel_phen", "B_gain_AFmax", "B_loss_AFmax", "DDD_disease", "OMIM_phenotype")
