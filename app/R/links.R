# External links of a gene and a variant, from Program_files/links.tsv (editable by the lab).

genome_build <- function() setting("genome_build", "PRIORR_GENOME_BUILD", "GRCh38")

read_links <- function(path = file.path(data_dir, "links.tsv")) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  links <- utils::read.delim(text = lines[!startsWith(lines, "#")], colClasses = "character", quote = "")
  links[nzchar(links$url), ]
}

# Values for the url templates. The variant is taken from variant_id (VCF representation) when present.
link_values <- function(v) {
  build <- genome_build()
  hg <- if (build == "GRCh37") "19" else "38"
  key <- if (!is.null(v$variant_id) && !is.na(v$variant_id)) strsplit(as.character(v$variant_id), "_", fixed = TRUE)[[1]]
         else as.character(c(v$CHROM, v$POS, v$REF, v$ALT))
  key <- c(key, rep(NA_character_, 4))[1:4]
  chrom <- key[1]
  if (!is.na(chrom) && !startsWith(chrom, "chr")) chrom <- paste0("chr", chrom)
  pos <- suppressWarnings(as.integer(key[2]))
  rsid <- regmatches(as.character(v$Existing_variation %||% ""), regexpr("rs[0-9]+", as.character(v$Existing_variation %||% "")))
  clean <- function(x) if (length(x) == 0 || is.na(x) || !nzchar(x)) NA_character_ else as.character(x)
  list(gene = clean(v$SYMBOL), chrom = clean(chrom), chrom_number = clean(sub("^chr", "", chrom)),
       pos = clean(pos), pos_from = clean(pos - 10), pos_to = clean(pos + 10),
       ref = clean(key[3]), alt = clean(key[4]), rsid = clean(rsid[1]), build = build,
       ucsc_build = paste0("hg", hg), ucsc_number = hg, gnomad_dataset = if (build == "GRCh37") "gnomad_r2_1" else "gnomad_r4")
}

# Values for the CNV url templates (AnnotSV columns)
cnv_link_values <- function(cnv) {
  values <- link_values(list())
  chrom <- sub("^chr", "", as.character(cnv$SV_chrom))
  c(values[c("build", "ucsc_build", "ucsc_number")],
    list(chrom = paste0("chr", chrom), chrom_number = chrom, start = format(as.numeric(cnv$SV_start), scientific = FALSE),
         end = format(as.numeric(cnv$SV_end), scientific = FALSE), build_lower = tolower(values$build),
         gnomad_sv_dataset = if (values$build == "GRCh37") "gnomad_sv_r2_1" else "gnomad_sv_r4"))
}

# Links of one type ("gene", "variant" or "cnv"): data.frame(name, url, description)
variant_links <- function(v, type, links = read_links()) {
  values <- if (type == "cnv") cnv_link_values(v) else link_values(v)
  links <- links[links$type == type, ]
  urls <- vapply(links$url, function(u) {
    for (k in names(values)) {
      if (!grepl(paste0("{", k, "}"), u, fixed = TRUE)) next
      if (is.na(values[[k]])) return(NA_character_)
      u <- gsub(paste0("{", k, "}"), utils::URLencode(values[[k]], reserved = TRUE), u, fixed = TRUE)
    }
    u
  }, "")
  data.frame(name = links$name, url = unname(urls), description = links$description, stringsAsFactors = FALSE)[!is.na(urls), ]
}
