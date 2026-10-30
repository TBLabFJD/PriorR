# Exon / intron structure of the transcript of a variant (Variant information > Gene structure).
# The exons of the transcripts (RefSeq NM_ / NR_ / XM_ / XR_ and Ensembl ENST of GENCODE) come from
# the tables of the UCSC Genome Browser, downloaded once into the shared folder for the whole lab
# (scripts/get_transcripts.R, started in the background the first time they are needed): no
# request per variant. Until then (or without internet) the EXON / INTRON columns of the file give a
# schematic drawing.

ucsc_tables <- function(build = genome_build()) {
  if (build == "GRCh37") list(genome = "hg19", tables = c("ncbiRefSeq", "wgEncodeGencodeCompV19"))
  else list(genome = "hg38", tables = c("ncbiRefSeq", "knownGene"))
}
transcript_db_path <- function(build = genome_build()) file.path(shared_dir, "annotations", paste0("transcripts_", build, ".parquet"))

ucsc_chrom <- function(chrom) {
  chrom <- sub("^chr", "", as.character(chrom))
  paste0("chr", if (chrom %in% c("M", "MT")) "M" else chrom)
}

# One transcript (a row of a genePred table) as list(name, gene, chrom, strand, tx_start, tx_end,
# cds_start, cds_end, exon_starts, exon_ends) with 0-based starts and 1-based ends, in genomic order
normalize_ucsc <- function(r) {
  nums <- function(x) as.numeric(strsplit(sub(",$", "", as.character(x)), ",", fixed = TRUE)[[1]])
  list(name = as.character(r$name), gene = as.character(r$name2 %||% ""), chrom = as.character(r$chrom),
       strand = as.character(r$strand), tx_start = as.numeric(r$txStart), tx_end = as.numeric(r$txEnd),
       cds_start = as.numeric(r$cdsStart), cds_end = as.numeric(r$cdsEnd),
       exon_starts = nums(r$exonStarts), exon_ends = nums(r$exonEnds))
}

# A genePred table of UCSC: ncbiRefSeq and GENCODE (16 columns: bin first, gene in name2) or
# knownGene (12 columns, without gene)
read_genepred <- function(path) {
  # DuckDB reads the .gz tables (fread would need the R.utils package)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  df <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_csv(%s, delim = '\t', header = false, all_varchar = true, quote = '', escape = '', null_padding = true)",
                                     sql_str(path)))
  if (ncol(df) >= 16) df <- df[, -1]
  data.frame(name = df[[1]], base = sub("\\.[0-9]+$", "", df[[1]]), chrom = df[[2]], strand = df[[3]],
             txStart = as.numeric(df[[4]]), txEnd = as.numeric(df[[5]]), cdsStart = as.numeric(df[[6]]),
             cdsEnd = as.numeric(df[[7]]), exonStarts = df[[9]], exonEnds = df[[10]],
             name2 = if (ncol(df) >= 15) df[[12]] else "", stringsAsFactors = FALSE)
}

# Downloads the transcript tables of UCSC for a genome build and writes the Parquet of the shared
# folder (scripts/get_transcripts.R)
build_transcript_db <- function(build = genome_build(), out = transcript_db_path(build)) {
  t <- ucsc_tables(build)
  old <- options(timeout = 900)
  on.exit(options(old))
  parts <- lapply(t$tables, function(table) {
    tmp <- tempfile(fileext = ".txt.gz")
    on.exit(unlink(tmp))
    utils::download.file(sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/database/%s.txt.gz", t$genome, table),
                         tmp, mode = "wb", quiet = TRUE)
    read_genepred(tmp)
  })
  write_transcript_db(do.call(rbind, parts), out)
}

write_transcript_db <- function(df, out) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  duckdb::duckdb_register(con, "t", df)
  tmp <- paste0(out, ".tmp", Sys.getpid())
  DBI::dbExecute(con, sprintf("COPY (SELECT * FROM t ORDER BY base) TO %s (FORMAT parquet, COMPRESSION zstd)", sql_str(tmp)))
  file.rename(tmp, out)
}

# Is the table of transcripts there? If not, it is downloaded in the background (once; the state
# in <table>.downloading / <table>.error). Returns "ready", "downloading" or "error: <message>".
transcript_db_state <- function(build = genome_build(), start = TRUE) {
  db <- transcript_db_path(build)
  if (file.exists(db)) return("ready")
  lock <- paste0(db, ".downloading")
  err <- paste0(db, ".error")
  if (file.exists(lock)) {
    l <- tryCatch(jsonlite::fromJSON(lock), error = function(e) list())
    if (process_alive(l$pid)) return("downloading")
    unlink(lock)
  }
  if (file.exists(err) && difftime(Sys.time(), file.info(err)$mtime, units = "mins") < 30)
    return(paste("error:", paste(readLines(err, warn = FALSE), collapse = " ")))
  if (!start) return("missing")
  dir.create(dirname(db), recursive = TRUE, showWarnings = FALSE)
  script <- file.path(data_dir, "scripts", "get_transcripts.R")
  system2(file.path(R.home("bin"), "Rscript"), c(shQuote(normalizePath(script)), build),
          stdout = paste0(db, ".log"), stderr = paste0(db, ".log"), wait = FALSE)
  "downloading"
}

# Why there is no drawing at scale (state of transcript_db_state)
transcript_db_note <- function(state, feature = NA) {
  if (identical(state, "downloading")) "the exon coordinates of the transcripts are being downloaded (once for the lab, ~30 MB, a minute or two)."
  else if (startsWith(state %||% "", "error")) paste("the exon coordinates could not be downloaded:", sub("^error: ", "", state))
  else sprintf("the transcript %s is not in the RefSeq / GENCODE tables of UCSC.", if (is.na(feature) || feature == "-") "" else feature)
}

# Rows of the table of transcripts (a SQL condition on name, base, name2, chrom)
transcript_rows <- function(where, db) {
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s) WHERE %s", sql_str(db), where))
}

# Exons of the transcript of the file (Feature): the same version, or another version of it, on the
# chromosome of the variant. NULL when it is not in the table (or the table is not there yet).
transcript_memo <- new.env()
transcript_structure <- function(feature, chrom, pos, build = genome_build(), db = transcript_db_path(build)) {
  if (is.null(feature) || is.na(feature) || !grepl("^(N[MR]_|X[MR]_|ENST)", feature)) return(NULL)
  key <- paste(db, feature, ucsc_chrom(chrom))
  if (exists(key, envir = transcript_memo, inherits = FALSE)) return(get(key, envir = transcript_memo))
  if (!file.exists(db)) return(NULL)
  rows <- transcript_rows(sprintf("base = %s AND chrom = %s", sql_str(sub("\\.[0-9]+$", "", feature)), sql_str(ucsc_chrom(chrom))), db)
  tx <- NULL
  if (nrow(rows)) {
    same <- rows[rows$name == feature, , drop = FALSE]
    if (!nrow(same)) same <- rows[order(-suppressWarnings(as.numeric(sub("^.*\\.", "", rows$name)))), , drop = FALSE]
    # the copy nearest to the variant (a few transcripts are in two places)
    same <- same[order(pmax(same$txStart - pos, pos - same$txEnd, 0)), , drop = FALSE]
    tx <- normalize_ucsc(as.list(same[1, ]))
  }
  assign(key, tx, envir = transcript_memo)
  tx
}

# Where the variant falls, in the order of the transcript (exon 1 = 5'):
# list(region = "exon" / "intron" / "upstream" / "downstream", number, total, size, from_start, to_end)
# from_start / to_end: bp from the 5' / 3' border of that exon or intron (1 = the first base)
variant_location <- function(tx, pos) {
  n <- length(tx$exon_starts)
  o <- if (tx$strand == "-") rev(seq_len(n)) else seq_len(n)
  s <- tx$exon_starts[o] + 1; e <- tx$exon_ends[o]   # 1-based, transcript order
  oriented <- function(a, b) if (tx$strand == "-") c(b - pos + 1, pos - a + 1) else c(pos - a + 1, b - pos + 1)
  i <- which(pos >= s & pos <= e)
  if (length(i)) {
    d <- oriented(s[i], e[i])
    return(list(region = "exon", number = i, total = n, size = e[i] - s[i] + 1, from_start = d[1], to_end = d[2]))
  }
  if (n > 1) for (j in seq_len(n - 1)) {
    lo <- if (tx$strand == "-") e[j + 1] + 1 else e[j] + 1
    hi <- if (tx$strand == "-") s[j] - 1 else s[j + 1] - 1
    if (pos >= lo && pos <= hi) {
      d <- oriented(lo, hi)
      return(list(region = "intron", number = j, total = n - 1, size = hi - lo + 1, from_start = d[1], to_end = d[2]))
    }
  }
  five <- if (tx$strand == "-") pos > max(e) else pos < min(s)
  list(region = if (five) "upstream" else "downstream", number = NA, total = n, size = NA,
       from_start = NA, to_end = NA)
}

# EXON / INTRON columns of VEP ("3/10")
parse_rank <- function(x) {
  if (is.null(x) || is.na(x) || !grepl("^[0-9]+(-[0-9]+)?/[0-9]+$", x)) return(NULL)
  p <- strsplit(x, "/", fixed = TRUE)[[1]]
  list(number = as.integer(sub("-.*", "", p[1])), total = as.integer(p[2]))
}

svg_text <- function(x, y, label, size = 12, anchor = "middle", color = "#37474f", weight = "normal")
  sprintf('<text x="%.1f" y="%.1f" font-size="%d" text-anchor="%s" fill="%s" font-weight="%s">%s</text>',
          x, y, size, anchor, color, weight, htmltools::htmlEscape(label))

# Drawing of the transcript in its own direction (5' left): exons at scale (thin UTR, thick coding
# part), introns shortened (square root of their length) so the whole gene fits; the variant in red
transcript_svg <- function(tx, pos, loc, width = 1000) {
  n <- length(tx$exon_starts)
  minus <- tx$strand == "-"
  o <- if (minus) rev(seq_len(n)) else seq_len(n)
  s <- tx$exon_starts[o]; e <- tx$exon_ends[o]   # 0-based half open, transcript order
  # segments in transcript order: exon 1, intron 1, exon 2...
  seg_from <- seg_to <- seg_w <- numeric(0); seg_exon <- logical(0)
  for (i in seq_len(n)) {
    seg_from <- c(seg_from, s[i]); seg_to <- c(seg_to, e[i]); seg_exon <- c(seg_exon, TRUE)
    seg_w <- c(seg_w, e[i] - s[i])
    if (i < n) {
      a <- if (minus) e[i + 1] else e[i]; b <- if (minus) s[i] else s[i + 1]
      seg_from <- c(seg_from, a); seg_to <- c(seg_to, b); seg_exon <- c(seg_exon, FALSE)
      seg_w <- c(seg_w, sqrt(max(b - a, 1)) * 4)
    }
  }
  left <- 40; right <- width - 40
  seg_w[seg_exon] <- pmax(seg_w[seg_exon], 3)
  scale <- (right - left) / sum(seg_w)
  x0 <- left + c(0, cumsum(seg_w * scale))[seq_along(seg_w)]
  # genomic position (0-based coordinate) -> x
  gx <- function(g) {
    k <- which(g >= seg_from & g <= seg_to)[1]
    if (is.na(k)) return(if ((g < min(seg_from)) != minus) left else right)
    f <- if (seg_to[k] > seg_from[k]) (g - seg_from[k]) / (seg_to[k] - seg_from[k]) else 0.5
    x0[k] + (if (minus) 1 - f else f) * seg_w[k] * scale
  }
  y <- 46
  parts <- c(sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#90a4ae" stroke-width="2"/>', left, y, right, y))
  for (i in seq_len(n)) {
    xa <- gx(s[i]); xb <- gx(e[i]); if (xb < xa) { t <- xa; xa <- xb; xb <- t }
    here <- identical(loc$region, "exon") && loc$number == i
    fill <- if (here) "#e65100" else "#26a69a"
    parts <- c(parts, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="12" fill="%s"><title>Exon %d (%d bp)</title></rect>',
                              xa, y - 6, max(xb - xa, 1.5), fill, i, as.integer(e[i] - s[i])))
    # coding part, thicker
    cs <- max(s[i], tx$cds_start); ce <- min(e[i], tx$cds_end)
    if (ce > cs) {
      ca <- gx(cs); cb <- gx(ce); if (cb < ca) { t <- ca; ca <- cb; cb <- t }
      parts <- c(parts, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="24" rx="2" fill="%s"><title>Exon %d (%d bp)</title></rect>',
                                ca, y - 12, max(cb - ca, 1.5), fill, i, as.integer(e[i] - s[i])))
    }
  }
  # exon numbers where there is room (always the first, the last and the one of the variant)
  last_x <- -Inf
  for (i in seq_len(n)) {
    xm <- (gx(s[i]) + gx(e[i])) / 2
    must <- i %in% c(1, n) || (identical(loc$region, "exon") && loc$number == i)
    if (must || xm - last_x > 22) {
      parts <- c(parts, svg_text(xm, y + 30, i, size = 11, color = if (must) "#37474f" else "#78909c"))
      last_x <- xm
    }
  }
  if (identical(loc$region, "intron")) {
    j <- loc$number
    a <- gx(if (minus) e[j + 1] else e[j]); b <- gx(if (minus) s[j] else s[j + 1])
    parts <- c(parts, sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e65100" stroke-width="4"/>', min(a, b), y, max(a, b), y))
  }
  xv <- gx(pos - 1)
  parts <- c(parts,
             sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#c62828" stroke-width="2"/>', xv, 10, xv, y + 14),
             sprintf('<polygon points="%.1f,10 %.1f,0 %.1f,0" fill="#c62828"/>', xv, xv - 6, xv + 6),
             svg_text(14, y + 5, "5′", size = 12, color = "#607d8b"),
             svg_text(width - 14, y + 5, "3′", size = 12, color = "#607d8b"))
  HTML(sprintf('<svg viewBox="0 0 %d 84" width="100%%" style="max-height: 110px; display: block;" role="img" aria-label="Exons of %s">%s</svg>',
               width, htmltools::htmlEscape(tx$name), paste(parts, collapse = "")))
}

cnv_color <- function(type) if (identical(toupper(type), "DUP")) "#1565c0" else "#c62828"

# Exons of a transcript (numbered from 5') inside a CNV (1-based start / end)
cnv_location <- function(tx, start, end) {
  n <- length(tx$exon_starts)
  o <- if (tx$strand == "-") rev(seq_len(n)) else seq_len(n)
  s <- tx$exon_starts[o]; e <- tx$exon_ends[o]
  exons <- which(e > start - 1 & s < end)
  list(exons = exons, total = n, whole = start - 1 <= tx$tx_start && end >= tx$tx_end)
}

# Transcript of each gene of a CNV: the one of the file (Tx) when given, otherwise the RefSeq NM_
# with more exons of the gene (on the chromosome of the CNV). NULL for a gene without data.
cnv_gene_transcripts <- function(genes, chrom, start, end, tx_names = list(), build = genome_build(),
                                 db = transcript_db_path(build)) {
  out <- stats::setNames(vector("list", length(genes)), genes)
  if (!length(genes) || !file.exists(db)) return(out)
  rows <- transcript_rows(sprintf("name2 IN (%s) AND chrom = %s", paste(sql_str(genes), collapse = ", "), sql_str(ucsc_chrom(chrom))), db)
  for (g in genes) {
    feature <- tx_names[[g]] %||% NA
    mine <- rows[rows$name2 == g, , drop = FALSE]
    if (!is.na(feature) && any(mine$base == sub("\\.[0-9]+$", "", feature))) mine <- mine[mine$base == sub("\\.[0-9]+$", "", feature), , drop = FALSE]
    else if (any(grepl("^N[MR]_", mine$name))) mine <- mine[grepl("^N[MR]_", mine$name), , drop = FALSE]
    if (!nrow(mine)) next
    n_exons <- lengths(strsplit(mine$exonStarts, ",", fixed = TRUE))
    out[g] <- list(normalize_ucsc(as.list(mine[which.max(n_exons), ])))
  }
  out
}

# One line for a CNV: the region of the CNV (and a margin) at genomic scale, the genes with their
# exons, the CNV shaded (red DEL, blue DUP) and the exons inside it in that color
cnv_region_svg <- function(txs, start, end, type, width = 1000) {
  txs <- Filter(Negate(is.null), txs)
  pad <- max((end - start) * 0.25, 2000)
  from <- start - pad; to <- end + pad
  left <- 40; right <- width - 40; y <- 46
  x <- function(g) left + (min(max(g, from), to) - from) / (to - from) * (right - left)
  color <- cnv_color(type)
  parts <- c(sprintf('<rect x="%.1f" y="16" width="%.1f" height="44" fill="%s" fill-opacity="0.13" stroke="%s" stroke-dasharray="4 3"/>',
                     x(start - 1), max(x(end) - x(start - 1), 2), color, color),
             sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cfd8dc" stroke-width="1"/>', left, y, right, y))
  last_label <- -Inf
  for (tx in txs[order(vapply(txs, function(t) t$tx_start, 1))]) {
    parts <- c(parts, sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#78909c" stroke-width="2"/>',
                              x(tx$tx_start), y, x(tx$tx_end), y))
    for (i in seq_along(tx$exon_starts)) {
      s <- tx$exon_starts[i]; e <- tx$exon_ends[i]
      if (e < from || s > to) next
      inside <- e > start - 1 && s < end
      parts <- c(parts, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="14" fill="%s"/>',
                                x(s), y - 7, max(x(e) - x(s), 1.5), if (inside) color else "#26a69a"))
    }
    xm <- (max(x(tx$tx_start), left) + min(x(tx$tx_end), right)) / 2
    if (xm - last_label > 60) {
      parts <- c(parts, svg_text(xm, y - 14, tx$gene, size = 11, color = "#37474f", weight = "bold"))
      last_label <- xm
    }
  }
  scale <- function(v) if (v >= 1e6) sprintf("%.1f Mb", v / 1e6) else sprintf("%.0f kb", v / 1e3)
  parts <- c(parts, svg_text(left, 78, format(round(from), big.mark = ","), size = 10, anchor = "start", color = "#90a4ae"),
             svg_text(right, 78, format(round(to), big.mark = ","), size = 10, anchor = "end", color = "#90a4ae"),
             svg_text((x(start - 1) + x(end)) / 2, 78, paste(toupper(type), scale(end - start + 1)), size = 10, color = color, weight = "bold"))
  HTML(sprintf('<svg viewBox="0 0 %d 84" width="100%%" style="max-height: 110px; display: block;" role="img" aria-label="CNV and genes">%s</svg>',
               width, paste(parts, collapse = "")))
}

# Short text of the genes of a CNV: "DMD: exons 38–45 of 79 · PMP22: whole gene"
cnv_genes_text <- function(txs, start, end, max_genes = 8) {
  txs <- Filter(Negate(is.null), txs)
  items <- vapply(txs, function(tx) {
    loc <- cnv_location(tx, start, end)
    what <- if (loc$whole) "whole gene" else if (!length(loc$exons)) "no exons" else {
      r <- range(loc$exons)
      sprintf("exon%s %s of %d", if (r[1] == r[2]) "" else "s", if (r[1] == r[2]) r[1] else paste0(r[1], "\u2013", r[2]), loc$total)
    }
    paste0(tx$gene, ": ", what)
  }, "")
  if (length(items) > max_genes) items <- c(items[seq_len(max_genes)], sprintf("+%d genes", length(items) - max_genes))
  paste(items, collapse = " \u00b7 ")
}

# Without coordinates: exons of the same size numbered, the exon or intron of the variant marked
schematic_svg <- function(exon, intron, width = 1000) {
  r <- parse_rank(exon) %||% parse_rank(intron)
  if (is.null(r)) return(NULL)
  in_exon <- !is.null(parse_rank(exon))
  n <- if (in_exon) r$total else r$total + 1
  left <- 40; right <- width - 40; y <- 40
  step <- (right - left) / n; w <- min(step * 0.55, 60)
  cx <- left + step * (seq_len(n) - 0.5)
  parts <- sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#90a4ae" stroke-width="2"/>', left, y, right, y)
  for (i in seq_len(n)) {
    here <- in_exon && i == r$number
    parts <- c(parts, sprintf('<rect x="%.1f" y="%d" width="%.1f" height="20" rx="2" fill="%s"/>', cx[i] - w / 2, y - 10, w,
                              if (here) "#e65100" else "#26a69a"))
    if (n <= 40 || i %in% c(1, n) || here || i %% 5 == 0) parts <- c(parts, svg_text(cx[i], y + 28, i, size = 11))
  }
  xv <- if (in_exon) cx[r$number] else (cx[r$number] + cx[r$number + 1]) / 2
  if (!in_exon) parts <- c(parts, sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e65100" stroke-width="4"/>',
                                          cx[r$number] + w / 2, y, cx[r$number + 1] - w / 2, y))
  parts <- c(parts,
             sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#c62828" stroke-width="2"/>', xv, 8, xv, y + 12),
             sprintf('<polygon points="%.1f,8 %.1f,0 %.1f,0" fill="#c62828"/>', xv, xv - 6, xv + 6),
             svg_text(14, y + 5, "5′", color = "#607d8b"), svg_text(width - 14, y + 5, "3′", color = "#607d8b"))
  HTML(sprintf('<svg viewBox="0 0 %d 76" width="100%%" style="max-height: 100px; display: block;" role="img" aria-label="Exons">%s</svg>',
               width, paste(parts, collapse = "")))
}

# Text of the location: "Exon 3 of 10 (156 bp): base 45 from its 5' end, 112 from its 3' end"
location_text <- function(loc) {
  n <- function(x) format(x, big.mark = ",")
  switch(loc$region,
    exon = sprintf("Exon %d of %d (%s bp): base %s from its 5′ end, %s from its 3′ end",
                   loc$number, loc$total, n(loc$size), n(loc$from_start), n(loc$to_end)),
    intron = sprintf("Intron %d of %d (%s bp): %s bp after exon %d, %s bp before exon %d",
                     loc$number, loc$total, n(loc$size), n(loc$from_start), loc$number, n(loc$to_end), loc$number + 1),
    upstream = "Before the first exon (5′ of the transcript)",
    downstream = "After the last exon (3′ of the transcript)")
}
