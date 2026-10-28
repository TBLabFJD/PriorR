ped_file <- function(lines) {
  f <- tempfile(fileext = ".ped")
  writeLines(lines, f)
  read_ped(f)
}

family <- ped_file(c("F P001 P002 P003 1 2", "F P002 0 0 1 1",
                     "F P003 0 0 2 1", "F P004 P002 P003 1 1"))

test_that("the sample genotype is the first _GT column", {
  expect_equal(sample_gt_column(c("X1_GT", "X1_DV_GT", "X1_DP")), "X1_GT")
})

test_that("PED roles are found from the file itself", {
  fam <- ped_samples(family, c("P001_GT", "P002_GT", "P003_GT", "P004_GT"))
  expect_equal(fam$affected, "P001_GT")
  expect_equal(fam$father, "P002_GT")
  expect_equal(fam$mother, "P003_GT")
  expect_equal(fam$relatives, "P004_GT")
  expect_equal(fam$unaffected_males, c("P002_GT", "P004_GT"))
})

test_that("compound heterozygous needs two different variants, not two transcripts", {
  df <- make_variants(4, SYMBOL = c("A", "A", "B", "B"), POS = c(1, 1, 2, 3))
  out <- inheritance_filter(df, "Recessive-Compound heterozygous")
  expect_equal(unique(out$SYMBOL), "B")
})

test_that("compound heterozygous with parents needs one variant from each", {
  df <- make_variants(5, SYMBOL = c("A", "A", "B", "B", "B"), POS = 1:5,
                      `P002_GT` = c("0/1", "0/1", "0/1", "0/0", "0/1"),
                      `P003_GT` = c("0/0", "0/0", "0/0", "0/1", "0/1"),
                      `P004_GT` = "0/0")
  out <- inheritance_filter(df, "Recessive-Compound heterozygous", family)
  expect_equal(out$POS, c(3, 4))   # gene A only paternal; variant 5 in both parents excluded
})

test_that("dominant, recessive and X-linked with a PED file", {
  df <- make_variants(5, CHROM = c("chr1", "chr1", "chr1", "chrX", "chrX"),
                      OMIM_phenotype = c("D, Autosomal dominant", "D, Autosomal dominant", "R, Autosomal recessive", "X", "X"),
                      `P001_GT` = c("0/1", "0/1", "1/1", "1/1", "1/1"),
                      `P002_GT` = c("0/0", "0/1", "0/1", "0/0", "1/1"),
                      `P003_GT` = c("0/0", "0/0", "0/1", "0/1", "0/1"),
                      `P004_GT` = c("0/0", "0/0", "0/1", "0/0", "0/0"))
  expect_equal(inheritance_filter(df, "Dominant", family)$POS, df$POS[1])
  expect_equal(inheritance_filter(df, "Recessive-Homozygous", family)$POS, df$POS[3])
  expect_equal(inheritance_filter(df, "X-linked", family)$POS, df$POS[4])   # unaffected father carries #5
})

test_that("inheritance without PED uses the sample genotype", {
  df <- make_variants(3, `P001_GT` = c("0/1", "1/1", "1/1"),
                      OMIM_phenotype = c("Autosomal dominant", "Autosomal recessive", "Autosomal dominant"))
  expect_equal(inheritance_filter(df, "Dominant")$POS, df$POS[1])
  expect_equal(inheritance_filter(df, "Recessive-Homozygous")$POS, df$POS[2])
  expect_equal(nrow(inheritance_filter(df, "")), 3)
})

test_that("ROH rows of the proband", {
  df <- make_variants(3, `P001_ROH` = c("True", "False", NA))
  expect_equal(roh_rows(df, "P001"), c(TRUE, FALSE, FALSE))
  expect_equal(roh_rows(df), c(TRUE, FALSE, FALSE))
})

test_that("de novo: proband heterozygous, both parents clean reference with enough depth", {
  df <- make_variants(7, CHROM = c("chr1", "chr1", "chr1", "chr1", "chr1", "chrX", "chr1"))
  df$`P_GT` <- c("0/1", "0/1", "0/1", "0/1", "0/1", "1", "0/1")
  df$`P_DP` <- c(30, 30, 30, 5, 30, 30, 30)
  df$`P_AD` <- c("15,15", "15,15", "15,15", "2,3", "15_15", "0,30", "27,3")
  df$`F_GT` <- c("0/0", "0/1", "0/0", "0/0", "./.", "1", "0/0")
  df$`F_DP` <- 30
  df$`F_AD` <- c("30,0", "15,15", "28,2", "30,0", "", "0,30", "30,0")
  df$`M_GT` <- "0/0"
  df$`M_DP` <- 30
  df$`M_AD` <- "30,0"
  df$`P001_GT` <- NULL; df$`P001_DP` <- NULL
  ped <- trio_ped("P", "F", "M", "1")
  out <- inheritance_filter(df, "De novo", ped)
  # 1 ok; 2 father carrier; 3 father 2 alt reads; 4 proband depth 5; 5 father no call;
  # 6 chrX male hemizygous, mother clean -> de novo; 7 proband VAF 0.1
  expect_equal(out$POS, df$POS[c(1, 6)])
  # female proband: chrX needs both parents clean (father hemizygous carrier here)
  expect_equal(inheritance_filter(df, "De novo", trio_ped("P", "F", "M", "2"))$POS, df$POS[1])
  # without parents there is no de novo
  expect_equal(nrow(inheritance_filter(df, "De novo", trio_ped("P", "", "", "1"))), 0)
  expect_equal(nrow(inheritance_filter(df, "De novo")), 0)
})

test_that("a trio chosen in the sidebar works like a PED file", {
  ped <- trio_ped("P", "F", "M", "2")
  expect_equal(ped$V2, c("P", "F", "M"))
  expect_equal(ped$V3[1], "F"); expect_equal(ped$V4[1], "M")
  expect_equal(ped$V6, c("2", "1", "1"))
  expect_equal(nrow(trio_ped("P", "", "M")), 2)
})
