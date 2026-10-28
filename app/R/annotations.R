# Annotations that PriorR adds when the variant file doesn't have them:
# - gene constraint of gnomAD v4.1 (LOEUF, pLI, missense Z; Program_files/Constraint, scripts/update_constraint.R)
# - AlphaMissense (Cheng et al., Science 2023; CC BY 4.0) for missense SNVs, from the Parquet file made by
#   scripts/get_alphamissense.R (setting alphamissense_file)

.annotation_cache <- new.env()

constraint_table <- function() {
  if (is.null(.annotation_cache$constraint)) {
    f <- file.path(data_dir, "Constraint", "gnomad_v4.1_constraint.tsv.gz")
    .annotation_cache$constraint <- if (file.exists(f)) utils::read.delim(gzfile(f), quote = "", stringsAsFactors = FALSE)
                                    else data.frame(SYMBOL = character(0))
  }
  .annotation_cache$constraint
}

add_constraint <- function(df) {
  if (!"SYMBOL" %in% names(df) || any(grepl("loeuf|oe_lof_upper|lof\\.oe_ci\\.upper", names(df), ignore.case = TRUE))) return(df)
  t <- constraint_table()
  if (!nrow(t)) return(df)
  i <- match(df$SYMBOL, t$SYMBOL)
  for (col in c("gnomAD_LOEUF", "gnomAD_pLI", "gnomAD_mis_z")) df[[col]] <- t[[col]][i]
  df
}

has_alphamissense <- function(df) any(grepl("alphamissense|am_pathogenicity", names(df), ignore.case = TRUE))

# AlphaMissense score and class of the SNVs of the table (NA for the others)
add_alphamissense <- function(df, file = alphamissense_file) {
  if (has_alphamissense(df) || !file.exists(file) || !all(c("CHROM", "POS", "REF", "ALT") %in% names(df))) return(df)
  snv <- nchar(as.character(df$REF)) == 1 & nchar(as.character(df$ALT)) == 1
  keys <- unique(data.frame(CHROM = sub("^(chr)?", "chr", as.character(df$CHROM[snv])),
                            POS = as.integer(df$POS[snv]), REF = as.character(df$REF[snv]),
                            ALT = as.character(df$ALT[snv]), stringsAsFactors = FALSE))
  df$AlphaMissense <- NA_real_
  df$AlphaMissense_class <- NA_character_
  if (!nrow(keys)) return(df)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  duckdb::duckdb_register(con, "k", keys)
  am <- DBI::dbGetQuery(con, sprintf(
    "SELECT a.CHROM, a.POS, a.REF, a.ALT, a.am_pathogenicity, a.am_class FROM read_parquet('%s') a
     JOIN k ON a.CHROM = k.CHROM AND a.POS = k.POS AND a.REF = k.REF AND a.ALT = k.ALT",
    gsub("'", "''", normalizePath(file))))
  key <- function(chrom, pos, ref, alt) paste(sub("^(chr)?", "chr", chrom), pos, ref, alt)
  i <- match(key(df$CHROM, df$POS, df$REF, df$ALT), key(am$CHROM, am$POS, am$REF, am$ALT))
  df$AlphaMissense <- round(am$am_pathogenicity[i], 3)
  df$AlphaMissense_class <- am$am_class[i]
  df
}

# All the annotations added by PriorR
add_annotations <- function(df) add_alphamissense(add_constraint(add_sgds(df)))
