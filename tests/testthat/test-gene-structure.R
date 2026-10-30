# Exon / intron structure of the transcript of a variant (gene_structure.R), synthetic transcripts

plus <- list(name = "NM_000001.2", gene = "GENE1", chrom = "chr1", strand = "+", tx_start = 100, tx_end = 600,
             cds_start = 150, cds_end = 550, exon_starts = c(100, 300, 500), exon_ends = c(200, 400, 600))
minus <- modifyList(plus, list(name = "ENST00000000001.3", strand = "-"))

# a small table of transcripts like the ones of UCSC: RefSeq (16 columns, gene in name2) and GENCODE
# knownGene (12 columns, without gene)
fake_transcript_db <- function() {
  dir <- tempfile(); dir.create(dir)
  refseq <- file.path(dir, "ncbiRefSeq.txt.gz"); known <- file.path(dir, "knownGene.txt.gz")
  rs <- c("585\tNM_000001.2\tchr1\t+\t100\t600\t150\t550\t3\t100,300,500,\t200,400,600,\t0\tGENE1\tcmpl\tcmpl\t0,0,0,",
          "585\tNM_000009.1\tchr1\t+\t100\t600\t150\t550\t1\t100,\t600,\t0\tGENE1\tcmpl\tcmpl\t0,",
          "585\tXM_000002.1\tchr1\t+\t2000\t2600\t2000\t2600\t2\t2000,2500,\t2100,2600,\t0\tGENE2\tcmpl\tcmpl\t0,0,",
          "585\tNM_000001.2\tchr1_KI270706v1_random\t+\t100\t600\t150\t550\t3\t100,300,500,\t200,400,600,\t0\tGENE1\tcmpl\tcmpl\t0,0,0,")
  writeLines(rs, gz <- gzfile(refseq)); close(gz)
  writeLines("ENST00000000001.3\tchr1\t-\t100\t600\t150\t550\t3\t100,300,500,\t200,400,600,\t\tuc001.1", gz <- gzfile(known)); close(gz)
  db <- file.path(dir, "transcripts_GRCh38.parquet")
  write_transcript_db(rbind(read_genepred(refseq), read_genepred(known)), db)
  db
}

test_that("the RefSeq and GENCODE tables of UCSC give the same structure", {
  db <- fake_transcript_db()
  rm(list = ls(transcript_memo), envir = transcript_memo)
  expect_equal(transcript_structure("NM_000001.2", "1", 150, db = db), plus)
  expect_equal(transcript_structure("ENST00000000001.3", "chr1", 150, db = db)[c("strand", "exon_starts")], list(strand = "-", exon_starts = c(100, 300, 500)))
  expect_equal(transcript_structure("ENST00000000001.3", "chr1", 150, db = db)$gene, "")
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

test_that("transcripts from the local table: another version when the file has a newer one, and nothing without the table", {
  db <- fake_transcript_db()
  rm(list = ls(transcript_memo), envir = transcript_memo)
  expect_equal(transcript_structure("NM_000001.3", "1", 150, db = db)$name, "NM_000001.2")
  expect_equal(transcript_structure("NM_000001.2", "1", 150, db = db)$chrom, "chr1")   # not the copy of the alt contig
  expect_null(transcript_structure("ENSR0001", "1", 150, db = db))
  expect_null(transcript_structure("NM_999.1", "1", 150, db = db))
  expect_null(transcript_structure("NM_000001.2", "1", 150, db = file.path(tempdir(), "none.parquet")))
  expect_equal(transcript_db_state(start = FALSE), "missing")   # shared folder of the tests: no table
  expect_match(transcript_db_note("downloading"), "being downloaded")
  rm(list = ls(transcript_memo), envir = transcript_memo)
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
  expect_equal(cnv_location(minus, 250, 450)$exons, 2L)
  expect_equal(cnv_location(minus, 150, 450)$exons, 2:3)   # minus strand: exon 3 is the first of the genome
  expect_true(cnv_location(plus, 50, 700)$whole)
  gene2 <- modifyList(plus, list(gene = "GENE2", tx_start = 2000, tx_end = 2600, exon_starts = c(2000, 2500), exon_ends = c(2100, 2600)))
  expect_equal(cnv_genes_text(list(plus, gene2, NULL), 250, 2650), "GENE1: exons 2\u20133 of 3 \u00b7 GENE2: whole gene")
  expect_equal(cnv_genes_text(list(plus), 210, 290), "GENE1: no exons")
  svg <- as.character(cnv_region_svg(list(plus, gene2), 250, 2650, "DUP"))
  expect_match(svg, "fill-opacity")
  expect_match(svg, "#1565c0")
  expect_match(svg, "GENE2")
  db <- fake_transcript_db()
  txs <- cnv_gene_transcripts(c("GENE1", "GENE2", "GENE3"), "1", 250, 450, db = db)
  expect_equal(txs$GENE1$name, "NM_000001.2")   # the NM_ with more exons
  expect_equal(txs$GENE2$name, "XM_000002.1")   # no NM_: another one
  expect_null(txs$GENE3)
  expect_equal(cnv_gene_transcripts("GENE1", "1", 250, 450, tx_names = list(GENE1 = "NM_000009"), db = db)$GENE1$name, "NM_000009.1")
})
