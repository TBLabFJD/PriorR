# Updates the Human Phenotype Ontology (terms, hierarchy and genes) in Program_files/HPO.
# Run from the app/ folder (needs internet):  Rscript Program_files/scripts/update_hpo.R
# or with the files already downloaded:       Rscript Program_files/scripts/update_hpo.R <folder>
# Source: https://hpo.jax.org (genes_to_phenotype.txt and hp.obo of the latest release)
args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
release <- "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
src <- if (length(args)) args[1] else {
  d <- tempfile("hpo-"); dir.create(d)
  for (f in c("genes_to_phenotype.txt", "hp.obo")) utils::download.file(file.path(release, f), file.path(d, f), mode = "wb", quiet = TRUE)
  d
}
out_dir <- file.path("Program_files", "HPO")
dir.create(out_dir, showWarnings = FALSE)

# Ontology: id, name and parents (is_a) of the terms that are not obsolete
obo <- readLines(file.path(src, "hp.obo"), encoding = "UTF-8")
version <- sub("^data-version: hp/releases/", "", grep("^data-version:", obo, value = TRUE)[1])
block <- cumsum(obo %in% c("[Term]", "[Typedef]"))
is_term <- block > 0 & block %in% block[obo == "[Term]"]
field <- function(prefix) {
  hit <- is_term & startsWith(obo, prefix)
  tapply(sub(paste0("^", prefix, " *"), "", obo[hit]), block[hit], function(x) x, simplify = FALSE)
}
ids <- unlist(field("id:")); names_ <- unlist(field("name:"))
parents <- field("is_a:"); obsolete <- unlist(field("is_obsolete:"))
terms <- data.frame(block = as.integer(names(ids)), id = ids, stringsAsFactors = FALSE)
terms$name <- names_[as.character(terms$block)]
terms$parents <- vapply(as.character(terms$block), function(b) paste(sub(" .*", "", parents[[b]] %||% character(0)), collapse = ";"), "")
terms <- terms[!as.character(terms$block) %in% names(obsolete)[obsolete == "true"], c("id", "name", "parents")]

# Genes annotated to each term (directly)
g2p <- utils::read.delim(file.path(src, "genes_to_phenotype.txt"), quote = "", stringsAsFactors = FALSE)
genes <- unique(data.frame(id = g2p$hpo_id, gene = g2p$gene_symbol, stringsAsFactors = FALSE))
genes <- genes[genes$id %in% terms$id, ]

# Number of genes of each term including its subterms (shown in PriorR)
children <- split(rep(terms$id, lengths(strsplit(terms$parents, ";"))), unlist(strsplit(terms$parents, ";")))
direct <- split(match(genes$gene, unique(genes$gene)), genes$id)
memo <- new.env()
all_genes <- function(id) {
  if (!is.null(memo[[id]])) return(memo[[id]])
  g <- unique(c(direct[[id]], unlist(lapply(children[[id]], all_genes))))
  memo[[id]] <- g
  g
}
terms$n_genes <- vapply(terms$id, function(id) length(all_genes(id)), 0L)

write_gz <- function(df, file) {
  con <- gzfile(file.path(out_dir, file), "w")
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}
write_gz(terms, "hpo_terms.tsv.gz")
write_gz(genes, "hpo_genes.tsv.gz")
writeLines(version, file.path(out_dir, "hpo_version.txt"))
message(sprintf("HPO %s: %d terms, %d with genes, %d genes", version, nrow(terms), sum(terms$n_genes > 0), length(unique(genes$gene))))
