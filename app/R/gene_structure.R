# Exon / intron structure of the transcript of a variant (Variant information > Gene structure).
# The exons of a transcript (NM_ / NR_ of RefSeq, ENST of Ensembl) come from the UCSC Genome Browser
# API the first time (only the name of the transcript and the position are sent) and are kept in the
# shared folder for everybody. Without them (no internet, unknown transcript) the EXON / INTRON
# columns of the file give a schematic drawing.

transcript_cache_dir <- function() file.path(shared_dir, "annotations", "transcripts")

ucsc_tracks <- function(build = genome_build()) {
  if (build == "GRCh37") list(genome = "hg19", tracks = c("ncbiRefSeq", "wgEncodeGencodeCompV19"))
  else list(genome = "hg38", tracks = c("ncbiRefSeq", "knownGene"))
}

ucsc_chrom <- function(chrom) {
  chrom <- sub("^chr", "", as.character(chrom))
  paste0("chr", if (chrom %in% c("M", "MT")) "M" else chrom)
}

# One transcript of the UCSC API (genePred of ncbiRefSeq / GENCODE, or bigGenePred of knownGene) as
# list(name, gene, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends)
# with 0-based starts and 1-based ends (BED-like), in genomic order
normalize_ucsc <- function(r) {
  nums <- function(x) as.numeric(strsplit(sub(",$", "", as.character(x)), ",", fixed = TRUE)[[1]])
  if (!is.null(r$exonStarts)) {
    starts <- nums(r$exonStarts); ends <- nums(r$exonEnds)
    tx <- c(r$txStart, r$txEnd); cds <- c(r$cdsStart, r$cdsEnd); gene <- r$name2
  } else {
    starts <- r$chromStart + nums(r$chromStarts); ends <- starts + nums(r$blockSizes)
    tx <- c(r$chromStart, r$chromEnd); cds <- c(r$thickStart, r$thickEnd); gene <- r$geneName
  }
  list(name = as.character(r$name), gene = as.character(gene %||% ""), chrom = as.character(r$chrom),
       strand = as.character(r$strand), tx_start = as.numeric(tx[1]), tx_end = as.numeric(tx[2]),
       cds_start = as.numeric(cds[1]), cds_end = as.numeric(cds[2]), exon_starts = starts, exon_ends = ends)
}

# Transcripts of UCSC that overlap a position (NULL without connection)
ucsc_transcripts <- function(chrom, pos, build = genome_build()) {
  t <- ucsc_tracks(build)
  old <- options(timeout = 8)
  on.exit(options(old))
  out <- list()
  for (track in t$tracks) {
    url <- sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=%s;chrom=%s;start=%d;end=%d",
                   t$genome, track, ucsc_chrom(chrom), as.integer(pos) - 1L, as.integer(pos))
    res <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(res)) return(if (length(out)) out else NULL)
    out <- c(out, lapply(res[[track]], normalize_ucsc))
  }
  out
}

# The transcript (Feature of the file) with the same version, or another version of it
pick_transcript <- function(records, feature) {
  base <- function(x) sub("\\.[0-9]+$", "", x)
  names <- vapply(records, `[[`, "", "name")
  i <- match(feature, names)
  if (is.na(i)) i <- match(base(feature), base(names))
  if (is.na(i)) NULL else records[[i]]
}

# Exons of a transcript: from the cache of the shared folder, or from UCSC (then cached). Not
# found or no connection: NULL (and, without connection, it isn't asked again in this session).
ucsc_offline <- new.env()
transcript_structure <- function(feature, chrom, pos, build = genome_build(), cache_dir = transcript_cache_dir(),
                                 fetch = ucsc_transcripts) {
  if (is.null(feature) || is.na(feature) || !grepl("^(N[MR]_|X[MR]_|ENST)", feature)) return(NULL)
  file <- file.path(cache_dir, build, paste0(gsub("[^A-Za-z0-9._-]", "_", feature), ".json"))
  if (file.exists(file)) {
    tx <- tryCatch(jsonlite::fromJSON(file, simplifyVector = TRUE), error = function(e) NULL)
    if (!is.null(tx) && length(tx$exon_starts)) return(tx)
  }
  if (isTRUE(ucsc_offline$down) && difftime(Sys.time(), ucsc_offline$when, units = "mins") < 10) return(NULL)
  records <- fetch(chrom, pos, build)
  if (is.null(records)) {
    ucsc_offline$down <- TRUE
    ucsc_offline$when <- Sys.time()
    return(NULL)
  }
  tx <- pick_transcript(records, feature)
  if (!is.null(tx)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(file, ".tmp", Sys.getpid())
    tryCatch({
      jsonlite::write_json(tx, tmp, auto_unbox = TRUE, digits = NA)
      file.rename(tmp, file)
    }, error = function(e) unlink(tmp))
  }
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
