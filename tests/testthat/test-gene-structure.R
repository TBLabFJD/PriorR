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
