test_that("without parameters nothing is filtered, rows without gene included", {
  df <- make_variants(5, SYMBOL = c("A", NA, "B", NA, "C"), Consequence = c("x", NA, "y", "x", "x"))
  out <- filter_snv(df)
  expect_equal(nrow(out), 5)
  expect_equal(out$POS, df$POS)
})

test_that("CADD is compared as a number, not as text", {
  df <- make_variants(3, CADD_PHRED = c(9.5, 15, NA))
  out <- filter_snv(df, snv_params(min_cadd = "15"))
  expect_equal(out$CADD_PHRED, c(15, NA))
})

test_that("frequencies keep missing values", {
  df <- make_variants(3, gnomADg_AF = c(0.5, 0.001, NA))
  expect_equal(nrow(filter_snv(df, snv_params(max_af = list(gnomADg_AF = 0.01)))), 2)
})

test_that("regions, consequences and types keep empty values when NA is selected", {
  df <- make_variants(3, Genomic_region = c("EXONIC", NA, "INTRONIC"))
  expect_equal(nrow(filter_snv(df, snv_params(region = c("EXONIC", "NA")))), 2)
  expect_equal(nrow(filter_snv(df, snv_params(region = "EXONIC"))), 1)
})

test_that("predicted RefSeq transcripts are excluded only with XM_/XR_", {
  df <- make_variants(3, Feature = c("NM_1", "XM_2", "ENST0000X"))
  expect_equal(filter_snv(df, snv_params(exclude_predicted = TRUE))$Feature, c("NM_1", "ENST0000X"))
})

test_that("several typed genes", {
  df <- make_variants(3, SYMBOL = c("A", "B", "C"))
  expect_equal(filter_snv(df, snv_params(genes = split_genes("A, C")))$SYMBOL, c("A", "C"))
})

test_that("hiding columns never changes the rows, any column can be hidden", {
  df <- make_variants(6, CLNSIG = c("Pathogenic", NA, "Benign", NA, "Pathogenic", NA), CADD_PHRED = c(30, 5, 20, 25, 1, 40))
  p <- snv_params(min_cadd = 10, exclude_clinvar_benign = TRUE, prioritize = TRUE)
  all_columns <- filter_snv(df, modifyList(p, list(columns = names(df))))
  few_columns <- filter_snv(df, modifyList(p, list(columns = "SYMBOL")))
  expect_equal(all_columns$POS, few_columns$POS)
  expect_equal(attr(few_columns, "display_columns"), c("SYMBOL", "Variants_in_gene"))
  expect_equal(attr(all_columns, "priority"), attr(few_columns, "priority"))
  # nothing selected: the variant key, so the table is never empty
  no_columns <- filter_snv(df, modifyList(p, list(columns = character(0))))
  expect_equal(attr(no_columns, "display_columns"), c("CHROM", "POS", "REF", "ALT", "Variants_in_gene"))
})

test_that("the same parameters always give the same table", {
  df <- make_variants(50, CADD_PHRED = rep(c(5, 30, NA, 22, 18), 10), CLNSIG = rep(c(NA, "Pathogenic"), 25),
                      SYMBOL = rep(c("A", "B", NA, "C", "D"), 10))
  p <- snv_params(min_cadd = 15, prioritize = TRUE, exclude_clinvar_benign = TRUE)
  expect_identical(filter_snv(df, p), filter_snv(df, p))
})

test_that("variants of another test are hidden or highlighted, matching indels by VCF representation", {
  df <- make_variants(3)
  df$variant_id[2] <- "chr1_1099_AT_A"
  ids <- variant_ids(c("1", "1"), c(1000, 1099), c("A", "AT"), c("G", "A"))
  hidden <- filter_snv(df, snv_params(other_test_ids = ids, hide_other_test = TRUE))
  expect_equal(hidden$POS, df$POS[3])
  shown <- filter_snv(df, snv_params(other_test_ids = ids))
  expect_equal(attr(shown, "not_in_other_test"), c(FALSE, FALSE, TRUE))
})

test_that("panels, excluded panels, phenotypes and BED", {
  df <- make_variants(4, SYMBOL = c("A", "B", "C", NA), POS = c(10, 20, 30, 40))
  expect_equal(filter_snv(df, snv_params(panel_genes = c("A", "C")))$SYMBOL, c("A", "C"))
  expect_equal(filter_snv(df, snv_params(excluded_genes = "A"))$POS, c(20, 30, 40))
  expect_equal(filter_snv(df, snv_params(phenotype_genes = "B"))$SYMBOL, "B")
  bed <- data.frame(chrom = "chr1", start = 15, end = 30)
  expect_equal(filter_snv(df, snv_params(bed = bed))$POS, c(20, 30))
})

test_that("classifications and comments are added", {
  df <- make_variants(2)
  ann <- data.frame(variant_key = variant_key("chr1", 1100, "A", "G"), classification = "Pathogenic", n_comments = 3L)
  out <- filter_snv(df, snv_params(annotations = ann))
  expect_equal(out$Classification, c(NA, "Pathogenic"))
  expect_equal(out$Comments, c(0L, 3L))
  expect_true(all(c("Classification", "Comments") %in% attr(out, "display_columns")))
})

test_that("CNV genes match whole names and sample filter works", {
  cnv <- data.frame(Samples_ID = c("S1", "S1", "S2"), Gene_name = c("BRCA12", "TP53,BRCA1", "BRCA1"),
                    SV_type = "DEL", ACMG_class = c("5", "3", "1"), N_PROGRAMS = 3, Annotation_mode = "full",
                    Location2 = c("exon", NA, "intron"), B_gain_AFmax = NA, B_loss_AFmax = c(0.5, NA, 0))
  out <- filter_cnv(cnv, cnv_params(samples = "S1", genes = "BRCA1"))
  expect_equal(out$Gene_name, "TP53,BRCA1")
  expect_equal(nrow(filter_cnv(cnv, cnv_params(min_acmg = 3, max_loss_af = 0.1))), 1)
  expect_equal(nrow(filter_cnv(cnv, cnv_params(region = c("exon", "NA")))), 2)
  # the rows that pass and those that don't, as a table like the filtered one
  kept <- filter_cnv(cnv, cnv_params(min_acmg = 3, columns = c("Gene_name", "ACMG_class")))
  expect_false(".source_row" %in% names(kept))
  expect_equal(attr(kept, "source_rows"), c(1L, 2L))
  ex <- excluded_cnv(cnv, 3L, cnv_params(columns = c("ACMG_class", "Gene_name")))
  expect_equal(ex$Gene_name, "BRCA1")
  expect_equal(attr(ex, "display_columns"), c("ACMG_class", "Gene_name"))
})

test_that("CNV files without some columns (AnnotSV without N_PROGRAMS) are not emptied", {
  cnv <- data.frame(Samples_ID = "S1", SV_type = c("DEL", "DUP", "DEL"), ACMG_class = c("3", "5", "1"),
                    AnnotSV_ranking_score = c(0.5, 1.9, -1), Annotation_mode = "full", Gene_name = "A")
  expect_equal(nrow(filter_cnv(cnv, cnv_params(min_programs = 3, min_acmg = 1))), 3)
  ordered <- filter_cnv(cnv, cnv_params(order_by_class = TRUE))
  expect_equal(ordered$ACMG_class, c("5", "3", "1"))
})

test_that("AnnotSV files: genotype of each sample and class of the split rows", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tSamples_ID\tS-001\tAnnotation_mode\tACMG_class",
               "1_100_200_DEL_1\t1\t100\t200\tDEL\tS-001\t0/1:PASS:357\tfull\t5",
               "1_100_200_DEL_1\t1\t100\t200\tDEL\tS-001\t0/1:PASS:357\tsplit\tfull=5"), f)
  cnv <- read_cnvs(f)
  expect_equal(cnv_samples(cnv), "S-001")
  expect_equal(cnv$`S-001_GT`, c("0/1", "0/1"))
  expect_equal(cnv$ACMG_class, c("5", "5"))
  links <- variant_links(cnv[1, ], "cnv")
  expect_equal(links$url[links$name == "UCSC"], "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr1:100-200")
})

test_that("HPO terms give the genes of the term and of its subterms", {
  expect_match(hpo_version(), "^20")
  all_genes <- hpo_genes("HP:0000001")
  rcd <- hpo_genes("HP:0000510")                         # Rod-cone dystrophy
  expect_true(all(c("RHO", "USH2A") %in% rcd))
  expect_true(all(rcd %in% hpo_genes("HP:0000478")))     # included in Abnormality of the eye
  expect_true(all(hpo_genes(c("HP:0000510", "HP:0001250")) %in% all_genes))
  expect_equal(unname(hpo_choices()[grepl("^Rod-cone dystrophy \\(HP:0000510\\)", names(hpo_choices()))]), "HP:0000510")
  expect_null(hpo_genes(character(0)))
})

test_that("HPO identifiers are read from pasted text", {
  expect_equal(parse_hpo_ids("HP:0000510, hp_0000556; 0001250 x"), c("HP:0000510", "HP:0000556", "HP:0001250"))
  expect_length(parse_hpo_ids("no ids"), 0)
})

test_that("frequencies or CADD that are not numbers count as missing values", {
  df <- make_variants(3, gnomADg_AF = c("", ".", "0.5"), CADD_PHRED = c("", "25", "."))
  expect_equal(nrow(filter_snv(df, snv_params(max_af = list(gnomADg_AF = 0.1)))), 2)
  expect_equal(nrow(filter_snv(df, snv_params(min_cadd = 20))), 3)
})

test_that("gnomAD constraint is added only when the file doesn't have LOEUF", {
  df <- add_constraint(make_variants(2, SYMBOL = c("SCN1A", "NOT_A_GENE")))
  expect_equal(df$gnomAD_LOEUF[1], 0.098)
  expect_true(is.na(df$gnomAD_LOEUF[2]))
  own <- make_variants(1, SYMBOL = "SCN1A", LOEUF = 0.5)
  expect_false("gnomAD_LOEUF" %in% names(add_constraint(own)))
})

test_that("AlphaMissense is added for SNVs only when the file doesn't have it", {
  f <- tempfile(fileext = ".parquet")
  con <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbExecute(con, sprintf("COPY (SELECT * FROM (VALUES ('chr1', 1000, 'A', 'G', 0.91::FLOAT, 'likely_pathogenic'),
                                     ('chr1', 1100, 'A', 'T', 0.05::FLOAT, 'likely_benign'))
                                     t(CHROM, POS, REF, ALT, am_pathogenicity, am_class)) TO '%s' (FORMAT parquet)", f))
  DBI::dbDisconnect(con, shutdown = TRUE)
  df <- make_variants(3, CHROM = c("chr1", "1", "chr1"), POS = c(1000, 1100, 1200), ALT = c("G", "T", "G"))
  out <- add_alphamissense(df, f)
  expect_equal(out$AlphaMissense, c(0.91, 0.05, NA))
  expect_equal(out$AlphaMissense_class[1], "likely_pathogenic")
  expect_identical(add_alphamissense(make_variants(1, AlphaMissense_score = 0.3), f)$AlphaMissense_score, 0.3)
  expect_false("AlphaMissense" %in% names(add_alphamissense(df, tempfile())))   # no file: nothing added
})

test_that("score filters: all of them (missing kept) or at least one", {
  df <- make_variants(4, REVELScore = c(0.9, 0.2, NA, NA), SpliceAI_DS_Max = c(0, 0.5, 0.6, NA))
  revel <- list(col = "REVELScore", op = ">=", value = 0.5)
  splice <- list(col = "SpliceAI_DS_Max", op = ">=", value = 0.2)
  expect_equal(nrow(filter_snv(df, snv_params(scores = list(revel)))), 3)            # 1, and 3-4 without REVEL
  expect_equal(filter_snv(df, snv_params(scores = list(revel, splice)))$POS, df$POS[3:4])  # missing scores are kept
  expect_equal(filter_snv(df, snv_params(scores = list(revel, splice), scores_mode = "any"))$POS, df$POS[1:3])
  loeuf <- list(col = "gnomAD_LOEUF", op = "<=", value = 0.6)
  expect_equal(nrow(filter_snv(cbind(df, gnomAD_LOEUF = c(0.2, 0.9, 0.5, NA)), snv_params(scores = list(loeuf)))), 3)
  expect_setequal(score_filters_for(c("REVELScore", "am_pathogenicity", "X"))$col, c("REVELScore", "am_pathogenicity"))
})

test_that("CNVs: GLOWgenes rank of a CNV is the best of its genes; filter and order by it", {
  df <- data.frame(SV_chrom = "1", SV_start = c(1, 100, 200), SV_end = c(50, 150, 250),
                   Gene_name = c("TNF;ABCA4", "TNF", "NOT_A_GENE"), stringsAsFactors = FALSE)
  src <- "orphanet:retinal_dystrophy_panel_glow"
  out <- filter_cnv(df, cnv_params(glow_sources = src, glow_order = TRUE))
  expect_equal(out$Gene_name[1], "TNF;ABCA4")              # ABCA4 is a panel gene: rank 0
  expect_equal(out$GLOW_rank[1], 0)
  expect_true(is.na(out$GLOW_rank[3]))
  expect_equal(attr(out, "display_columns")[1], "GLOW_rank")
  expect_equal(nrow(filter_cnv(df, cnv_params(glow_sources = src, glow_max = 0))), 1)
  expect_equal(cnv_locus(df[1, ]), "1:1-550")
  big <- data.frame(SV_chrom = "chr2", SV_start = 1e6, SV_end = 2e6)
  expect_equal(cnv_locus(big), "chr2:999700-1000300 chr2:1999700-2000300")
})

test_that("filtered out: the rows that don't pass, as a table like the filtered one", {
  df <- make_variants(6, CADD_PHRED = c(30, 5, 20, NA, 16, 40), CLNSIG = c(NA, NA, "Benign", NA, "Pathogenic", NA))
  p <- snv_params(min_cadd = 15, exclude_clinvar_benign = TRUE, prioritize = TRUE, columns = c("SYMBOL", "CADD_PHRED"))
  kept <- filter_snv(df, p)
  expect_false(".source_row" %in% names(kept))
  expect_setequal(attr(kept, "source_rows"), c(1, 4, 5, 6))
  out <- excluded_rows(df, kept)
  expect_equal(out, c(2, 3))
  ex <- excluded_snv(df, out, p)
  expect_equal(ex$CADD_PHRED, c(5, 20))
  expect_equal(attr(ex, "display_columns"), c("SYMBOL", "Variants_in_gene", "CADD_PHRED"))
  expect_equal(length(attr(ex, "priority")), 2)
  # pass + filtered out = all the rows
  expect_equal(nrow(kept) + length(out), nrow(df))
  # one row per variant (GLOWgenes view): another transcript of a shown variant is not out
  two <- df[c(1, 1, 2), ]
  kept1 <- filter_snv(two, snv_params(min_cadd = 15))
  kept1 <- kept1[1, ]
  attr(kept1, "source_rows") <- 1L
  expect_equal(excluded_rows(two, kept1), c(2, 3))
  expect_equal(excluded_rows(two, kept1, by_variant = TRUE), 3)
})

test_that("Variants_in_gene: different variants of the gene of each row that pass the filters", {
  df <- make_variants(5, SYMBOL = c("A", "A", "B", "A", NA), CADD_PHRED = c(30, 25, 20, 2, 40))
  df <- rbind(df, df[1, ])   # another transcript of the first variant
  kept <- filter_snv(df, snv_params(min_cadd = 10))
  expect_equal(kept$Variants_in_gene, c(2L, 2L, 1L, NA, 2L))
  # the filtered out rows count the variants of their gene that pass
  ex <- excluded_snv(df, excluded_rows(df, kept), kept = kept)
  expect_equal(ex$SYMBOL, "A")
  expect_equal(ex$Variants_in_gene, 2L)
})
