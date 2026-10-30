test_that("variant ids ignore the chr prefix and case", {
  expect_equal(variant_ids("chr1", 100, "a", "g"), variant_ids("1", "100", "A", "G"))
})

test_that("row ids use the VCF representation of variant_id when present", {
  df <- data.frame(CHROM = "chr10", POS = 13274236, REF = "T", ALT = NA, variant_id = "chr10_13274235_CT_C")
  expect_equal(row_variant_ids(df), variant_ids("10", 13274235, "CT", "C"))
  expect_equal(row_variant_ids(df[, 1:4]), variant_ids("chr10", 13274236, "T", NA))
})

test_that("empty values are kept when NA is selected", {
  expect_equal(in_selection(c("EXONIC", NA, "INTRONIC"), c("EXONIC", "NA")), c(TRUE, TRUE, FALSE))
  expect_equal(in_selection(c("EXONIC", NA), "EXONIC"), c(TRUE, FALSE))
})

test_that("typed genes are split", {
  expect_equal(split_genes(" BRCA1, BRCA2;TP53  ABCA4 "), c("BRCA1", "BRCA2", "TP53", "ABCA4"))
  expect_equal(split_genes(""), character(0))
})

test_that("gene fields with several genes match whole names only", {
  expect_equal(genes_overlap(c("BRCA12", "TP53,BRCA1", "A;B"), c("BRCA1", "B")), c(FALSE, TRUE, TRUE))
})

test_that("BED regions use 0-based starts and handle overlaps", {
  bed <- data.frame(chrom = c("chr1", "chr1", "2"), start = c(100, 150, 10), end = c(200, 160, 20))
  df <- data.frame(CHROM = c("1", "chr1", "chr1", "chr1", "chr2", "chr3"), POS = c(100, 101, 200, 201, 15, 15))
  expect_equal(in_bed(df, bed), c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE))
  # a long region followed by a short one inside it
  bed2 <- data.frame(chrom = "chr1", start = c(100, 120), end = c(1000, 130))
  expect_true(in_bed(data.frame(CHROM = "chr1", POS = 500), bed2))
})

test_that("external links of a gene and a variant (VCF representation; no rsID, no ClinVar link)", {
  v <- list(SYMBOL = "ERCC6", CHROM = "chr10", POS = 49472376, REF = "C", ALT = "T",
            Existing_variation = "rs145720191,COSV1", variant_id = "chr10_49472376_C_T")
  gene <- variant_links(v, "gene")
  expect_true(all(c("GeneCards", "STRING") %in% gene$name))
  expect_equal(gene$url[gene$name == "GeneCards"], "https://www.genecards.org/cgi-bin/carddisp.pl?gene=ERCC6")
  var <- variant_links(v, "variant")
  expect_equal(var$url[var$name == "gnomAD"], "https://gnomad.broadinstitute.org/variant/10-49472376-C-T?dataset=gnomad_r4")
  expect_equal(var$url[var$name == "dbSNP"], "https://www.ncbi.nlm.nih.gov/snp/rs145720191")
  deletion <- list(SYMBOL = "TYSND1", CHROM = "chr10", POS = 70146046, REF = "CTC", ALT = NA,
                   Existing_variation = NA, variant_id = "chr10_70146045_ACTC_A")
  var <- variant_links(deletion, "variant")
  expect_false("dbSNP" %in% var$name)
  expect_match(var$url[var$name == "VarSome"], "chr10-70146045-ACTC-A$")
})

test_that("reads: index next to BAM / CRAM, sample from the file name, range requests", {
  dir <- tempfile(); dir.create(dir)
  bam <- file.path(dir, "run1_P001.sorted.bam"); writeBin(as.raw(0:255), bam)
  expect_true(is.na(alignment_index(bam)))
  file.create(paste0(bam, ".bai"))
  expect_equal(alignment_index(bam), paste0(bam, ".bai"))
  cram <- file.path(dir, "x.cram"); file.create(cram, file.path(dir, "x.crai"))
  expect_equal(alignment_index(cram), file.path(dir, "x.crai"))
  expect_equal(file_sample(bam, c("P001", "10-04")), "P001")
  expect_equal(file_sample(cram, "P001"), "x")

  piece <- range_response(list(path = bam), list(HTTP_RANGE = "bytes=10-19"))
  expect_equal(piece$status, 206L)
  expect_equal(piece$content, as.raw(10:19))
  expect_equal(piece$headers$`Content-Range`, "bytes 10-19/256")
  open_end <- range_response(list(path = bam), list(HTTP_RANGE = "bytes=250-"))
  expect_equal(open_end$content, as.raw(250:255))
  expect_equal(range_response(list(path = bam), list())$status, 200L)
  expect_equal(range_response(list(path = bam), list(HTTP_RANGE = "bytes=300-400"))$status, 416L)
})

test_that("column labels: SYMBOL as GENE, no underscores, sample field first", {
  expect_equal(column_labels(c("SYMBOL", "gnomADg_AF_popmax", "P001_GT", "P001_DV_GT", "CHROM")),
               c("GENE", "gnomADg AF popmax", "GT P001", "GT P001 DV", "CHROM"))
  expect_equal(unname(column_choices(c("SYMBOL", "HGVSc"))), c("SYMBOL", "HGVSc"))
  expect_equal(names(column_choices(c("SYMBOL", "HGVSc"))), c("GENE", "HGVSc"))
})

test_that("column order: the picker lists the shown columns first, dragging a header moves them", {
  all <- c("CHROM", "POS", "REF", "ALT", "SYMBOL", "HGVSc", "CADD_PHRED")
  expect_equal(unname(ordered_column_choices(all, c("SYMBOL", "CHROM", "nope"))),
               c("SYMBOL", "CHROM", "POS", "REF", "ALT", "HGVSc", "CADD_PHRED"))
  # the table shows SYMBOL, CHROM, POS (+ GLOW_rank and Classify); CADD_PHRED is selected but not shown
  expect_equal(reordered_columns(c("SYMBOL", "CHROM", "POS", "CADD_PHRED"), c("GLOW_rank", "POS", "SYMBOL", "CHROM", "Classify")),
               c("POS", "SYMBOL", "CHROM", "CADD_PHRED"))
})

test_that("BAM / CRAM indexes are found with any naming, also in capitals", {
  dir <- tempfile(); dir.create(dir)
  f <- function(...) { p <- file.path(dir, ...); file.create(p); p }
  a <- f("P001.bam"); ai <- f("P001.bam.bai")
  b <- f("P002.BAM"); bi <- f("P002.BAI")
  c <- f("P003.cram"); ci <- f("P003.crai")
  d <- f("P004.bam")
  # macOS and Windows don't tell capitals from small letters in file names: compare ignoring case
  same_file <- function(x, y) expect_equal(tolower(normalizePath(x)), tolower(normalizePath(y)))
  same_file(alignment_index(a), ai)
  same_file(alignment_index(b), bi)
  same_file(alignment_index(c), ci)
  expect_true(is.na(alignment_index(d)))
  expect_equal(is_alignment_file(c("x.BAM", "x.cram", "x.bam.bai", "x.txt")), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("VCF of another test, plain or .gz (no zcat), with multiallelic sites split", {
  dir <- tempfile(); dir.create(dir)
  lines <- c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
             "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1", "chr2\t200\trs1\tC\tT,G\t50\tPASS\t.\tGT\t1/2")
  writeLines(lines, file.path(dir, "a.vcf"))
  con <- gzfile(file.path(dir, "a.vcf.gz"), "w"); writeLines(lines, con); close(con)
  ids <- read_vcf_ids(file.path(dir, "a.vcf"))
  expect_equal(ids, variant_ids(c("chr1", "chr2", "chr2"), c("100", "200", "200"), c("A", "C", "C"), c("G", "T", "G")))
  expect_identical(read_vcf_ids(file.path(dir, "a.vcf.gz")), ids)
})

test_that("patient in the header: proband, other samples, file as tooltip", {
  h <- as.character(patient_header("10-0432", c("P002", "P003"), "10-0432.tsv"))
  expect_match(h, "<b>10-0432</b>")
  expect_match(h, "\\+ P002, P003")
  expect_match(h, 'title="File: 10-0432.tsv"')
  expect_match(as.character(patient_header(character(0), file = "x.tsv")), "x.tsv")
  expect_null(patient_header(NA))
})

test_that("table order: empty frequencies sort as 0, earlier sorts break the ties", {
  expect_equal(is_af_column(c("gnomADg_AF", "gnomADg_AF_popmax", "FJD_MAF_AF", "CSVS_AF", "P001_VAF", "CADD_PHRED")),
               c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(sort_key(c(0.2, NA, 0.01), "gnomADe_AF"), c(0.2, 0, 0.01))
  expect_equal(sort_key(c(20, NA, 5), "CADD_PHRED"), c(20, NA, 5))
  expect_equal(sort_key(c("10", "9", NA), "x"), c(10, 9, NA))   # numbers as text
  df <- data.frame(SYMBOL = c("B", "A", "B", "A"), gnomADg_AF = c(0.1, NA, 0.001, 0.1), stringsAsFactors = FALSE)
  columns <- stats::setNames(lapply(names(df), function(n) list(name = n, orderable = "true", searchable = "true",
                                                                search = list(value = "", regex = "false"))), c("0", "1"))
  params <- function(order) list(draw = "1", start = "0", length = "-1", escape = "false", columns = columns,
                                 search = list(value = "", regex = "false", caseInsensitive = "true"), order = order)
  # gnomAD ascending: the empty one first
  res <- priorr_dt_filter(df, params(list(list(column = "1", dir = "asc"))))
  expect_equal(res$DT_rows_all, c(2L, 3L, 1L, 4L))
  # gene, then (clicked before) frequency descending within each gene
  res <- priorr_dt_filter(df, params(list(list(column = "0", dir = "asc"), list(column = "1", dir = "desc"))))
  expect_equal(res$DT_rows_all, c(4L, 2L, 1L, 3L))
})

test_that("a table drawn again keeps the filter boxes, search, order and page length by column name", {
  state <- table_view_state(c("CHROM", "SYMBOL", "CADD_PHRED"), c("", "LZ", "10 ... 40"), "chr1",
                            list(order = list(list("CADD_PHRED", "desc"), list("SYMBOL", "asc")), length = 25))
  expect_equal(state$search_cols, c(SYMBOL = "LZ", CADD_PHRED = "10 ... 40"))
  expect_null(table_view_state(NULL, "x", "", list()))
  # the new table has no CHROM and SYMBOL in another place
  df <- data.frame(SYMBOL = "A", POS = 1, CADD_PHRED = 2)
  o <- priorr_datatable(df, "full", state = state)$x$options
  expect_equal(o$searchCols, list(list(search = "LZ"), NULL, list(search = "10 ... 40")))
  expect_equal(o$search, list(search = "chr1"))
  expect_equal(o$order, list(list(2L, "desc"), list(0L, "asc")))
  expect_equal(o$pageLength, 25L)
  # nothing kept: DT defaults
  expect_null(priorr_datatable(df, "full")$x$options$searchCols)
})
