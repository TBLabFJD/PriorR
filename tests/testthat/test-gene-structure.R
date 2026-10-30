# Exon / intron structure of the transcript of a variant (gene_structure.R), synthetic transcripts

plus <- list(name = "NM_000001.2", gene = "GENE1", chrom = "chr1", strand = "+", tx_start = 100, tx_end = 600,
             cds_start = 150, cds_end = 550, exon_starts = c(100, 300, 500), exon_ends = c(200, 400, 600))
minus <- modifyList(plus, list(name = "ENST00000000001.3", strand = "-"))

test_that("UCSC records of RefSeq (genePred) and GENCODE (bigGenePred) give the same structure", {
  refseq <- normalize_ucsc(list(name = "NM_000001.2", name2 = "GENE1", chrom = "chr1", strand = "+", txStart = 100, txEnd = 600,
                                cdsStart = 150, cdsEnd = 550, exonStarts = "100,300,500,", exonEnds = "200,400,600,"))
  expect_equal(refseq, plus)
  gencode <- normalize_ucsc(list(name = "NM_000001.2", geneName = "GENE1", chrom = "chr1", strand = "+", chromStart = 100,
                                 chromEnd = 600, thickStart = 150, thickEnd = 550, blockSizes = "100,100,100,", chromStarts = "0,200,400,"))
  expect_equal(gencode, plus)
})

test_that("location of the variant in the order of the transcript, with distances to the borders", {
  expect_equal(variant_location(plus, 150)[c("region", "number", "total", "size", "from_start", "to_end")],
               list(region = "exon", number = 1L, total = 3L, size = 100, from_start = 50, to_end = 51))
  expect_equal(variant_location(plus, 250)[c("region", "number", "size", "from_start", "to_end")],
               list(region = "intron", number = 1L, size = 100, from_start = 50, to_end = 51))
  # minus strand: exon 1 is the last one of the genome, distances from its 5' end (the higher position)
  expect_equal(variant_location(minus, 550)[c("region", "number", "from_start", "to_end")],
               list(region = "exon", number = 1L, from_start = 51, to_end = 50))
  expect_equal(variant_location(minus, 450)[c("region", "number", "from_start", "to_end")],
               list(region = "intron", number = 1L, from_start = 51, to_end = 50))
  expect_equal(variant_location(minus, 50)$region, "downstream")
  expect_equal(variant_location(plus, 50)$region, "upstream")
  expect_match(location_text(variant_location(plus, 250)), "Intron 1 of 2 \\(100 bp\\): 50 bp after exon 1, 51 bp before exon 2")
})

test_that("exons are asked once and kept in the cache; another version of the transcript is used", {
  dir <- tempfile()
  calls <- 0
  fetch <- function(chrom, pos, build) { calls <<- calls + 1; list(modifyList(plus, list(name = "NM_000002.1")), plus) }
  rm(list = ls(ucsc_offline), envir = ucsc_offline)
  tx <- transcript_structure("NM_000001.2", "1", 150, cache_dir = dir, fetch = fetch)
  expect_equal(tx$name, "NM_000001.2")
  tx <- transcript_structure("NM_000001.2", "1", 150, cache_dir = dir, fetch = fetch)
  expect_equal(calls, 1)
  expect_equal(tx$exon_ends, c(200, 400, 600))
  expect_equal(transcript_structure("NM_000001.3", "1", 150, cache_dir = dir, fetch = fetch)$name, "NM_000001.2")
  # not a transcript, or no connection: nothing (and no more questions for a while)
  expect_null(transcript_structure("ENSR0001", "1", 150, cache_dir = dir, fetch = fetch))
  expect_null(transcript_structure("NM_999.1", "1", 150, cache_dir = dir, fetch = function(...) NULL))
  expect_null(transcript_structure("NM_998.1", "1", 150, cache_dir = dir, fetch = function(...) stop("asked again")))
  rm(list = ls(ucsc_offline), envir = ucsc_offline)
})

test_that("drawings: at scale with the variant, or schematic from EXON / INTRON", {
  svg <- as.character(transcript_svg(minus, 450, variant_location(minus, 450)))
  expect_match(svg, "<svg")
  expect_match(svg, "Exon 3 \\(100 bp\\)")
  expect_equal(parse_rank("3/10"), list(number = 3L, total = 10L))
  expect_null(parse_rank(NA))
  expect_match(as.character(schematic_svg("3/10", NA)), "<svg")
  expect_match(as.character(schematic_svg("-", "2/9")), "<svg")
  expect_null(schematic_svg("-", "-"))
})

test_that("CNVs: exons of the gene inside the CNV, and the transcript of each gene", {
  # plus strand exons 101-200, 301-400, 501-600: a deletion of 250-450 takes exon 2
  loc <- cnv_location(plus, 250, 450)
  expect_equal(loc$exons, 2L)
  expect_false(loc$whole)
  expect_equal(cnv_location_text(loc, "DEL"), "Deleted: exon 2 of 3 (part of the gene)")
  expect_equal(cnv_location(minus, 250, 450)$exons, 2L)
  expect_equal(cnv_location(minus, 150, 450)$exons, 2:3)   # minus strand: exon 3 is the first of the genome
  expect_match(cnv_location_text(cnv_location(plus, 50, 700), "DUP"), "Duplicated: the whole gene \\(3 exons\\)")
  expect_match(cnv_location_text(cnv_location(plus, 210, 290), "DEL"), "Inside an intron")
  svg <- as.character(transcript_svg(plus, 250, loc, cnv = list(start = 250, end = 450, type = "DEL")))
  expect_match(svg, "fill-opacity")
  expect_match(svg, "#c62828")

  dir <- tempfile()
  rm(list = ls(ucsc_offline), envir = ucsc_offline)
  calls <- 0
  other <- modifyList(plus, list(name = "NM_000009.1", gene = "GENE1", exon_starts = 100, exon_ends = 600))
  fetch <- function(chrom, pos, build, end = pos) { calls <<- calls + 1; list(other, plus, modifyList(plus, list(name = "XM_1", gene = "GENE2"))) }
  txs <- cnv_gene_transcripts(c("GENE1", "GENE2", "GENE3"), "1", 250, 450, cache_dir = dir, fetch = fetch)
  expect_equal(txs$GENE1$name, "NM_000001.2")   # the NM_ with more exons
  expect_equal(txs$GENE2$name, "XM_1")          # no NM_: another one
  expect_null(txs$GENE3)
  txs <- cnv_gene_transcripts(c("GENE1", "GENE2"), "1", 250, 450, cache_dir = dir, fetch = fetch)
  expect_equal(calls, 1)                         # kept by gene
  expect_equal(txs$GENE1$exon_ends, c(200, 400, 600))
  rm(list = ls(ucsc_offline), envir = ucsc_offline)
})
