# Human Phenotype Ontology (Program_files/HPO, made by scripts/update_hpo.R): the genes of a term
# are the genes annotated to it or to any of its subterms, so several terms work like a gene panel.

hpo_dir <- function() file.path(data_dir, "HPO")
.hpo_cache <- new.env()

hpo_data <- function() {
  if (is.null(.hpo_cache$terms)) {
    terms <- utils::read.delim(gzfile(file.path(hpo_dir(), "hpo_terms.tsv.gz")), quote = "", colClasses = "character")
    terms$n_genes <- as.integer(terms$n_genes)
    genes <- utils::read.delim(gzfile(file.path(hpo_dir(), "hpo_genes.tsv.gz")), quote = "", colClasses = "character")
    parents <- strsplit(terms$parents, ";")
    .hpo_cache$children <- split(rep(terms$id, lengths(parents)), unlist(parents))
    .hpo_cache$genes <- split(genes$gene, genes$id)
    .hpo_cache$terms <- terms
  }
  .hpo_cache
}

hpo_version <- function() {
  f <- file.path(hpo_dir(), "hpo_version.txt")
  if (file.exists(f)) readLines(f, n = 1) else NA_character_
}

# Choices of the selectize: "Name (HP:...) [n genes]" -> id, only terms with genes
hpo_choices <- function() {
  t <- hpo_data()$terms
  t <- t[t$n_genes > 0, ]
  stats::setNames(t$id, sprintf("%s (%s) [%d genes]", t$name, t$id, t$n_genes))
}

# Genes of the terms and of all their subterms
hpo_genes <- function(ids) {
  if (!length(ids)) return(NULL)
  h <- hpo_data()
  seen <- character(0)
  todo <- ids
  while (length(todo)) {
    seen <- c(seen, todo)
    todo <- setdiff(unlist(h$children[todo], use.names = FALSE), seen)
  }
  sort(unique(unlist(h$genes[seen], use.names = FALSE)))
}

# HPO identifiers in a text ("HP:0000510, hp_0000556 0001250..." -> "HP:0000510" ...)
parse_hpo_ids <- function(text) {
  ids <- regmatches(text, gregexpr("(?i)hp[:_]?[0-9]{7}|\\b[0-9]{7}\\b", text, perl = TRUE))[[1]]
  unique(sprintf("HP:%s", substring(ids, nchar(ids) - 6)))
}
