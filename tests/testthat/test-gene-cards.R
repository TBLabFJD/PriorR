# Cards of Variant information about the gene (gene_cards.R)

test_that("panels of the library with the gene and their PanelApp level", {
  catalog <- data.frame(id = c("a", "b", "c", "d"), source = c("PanelApp England", "PanelApp England", "Orphanet", "Saved by the users"),
                        name = c("Retina", "Deafness", "Usher", "Mine"), stringsAsFactors = FALSE)
  catalog$genes <- list("USH2A", c("MYO7A"), c("USH2A", "MYO7A"), character(0))
  catalog$amber <- list(character(0), "USH2A", character(0), character(0))
  catalog$red <- list(character(0), character(0), character(0), character(0))
  p <- gene_panels("USH2A", catalog)
  expect_equal(p$name, c("Retina", "Usher", "Deafness"))
  expect_equal(p$level, c("green", "in", "amber"))
  expect_equal(nrow(gene_panels("TP53", catalog)), 0)
  expect_equal(nrow(gene_panels(NA, catalog)), 0)
})

test_that("other variants of the gene, with their parental origin in a trio", {
  ped <- trio_ped("P001", "P002", "P003")
  df <- make_variants(4, SYMBOL = c("ABCA4", "ABCA4", "ABCA4", "RHO"), POS = c(100, 200, 200, 300),
                      HGVSp = c("NP_1:p.Arg1Ter", NA, NA, "NP_2:p.Pro23His"), HGVSc = c("NM_1:c.1C>T", "NM_1:c.2-5A>G", "NM_1:c.2-5A>G", "NM_2:c.3G>A"))
  df$P002_GT <- c("0/1", "0/0", "0/0", "0/1"); df$P003_GT <- c("0/0", "0/1", "0/1", "0/0")
  others <- gene_other_variants(df, df[1, ], ped)
  expect_equal(others$position, "chr1:200 A>G")       # the same variant twice (two transcripts) is one
  expect_equal(others$change, "c.2-5A>G")             # no HGVSp: HGVSc
  expect_equal(others$origin, "mother")
  expect_match(others$genotypes, "P001 0/1 · P002 0/0 · P003 0/1")
  expect_equal(variant_origin(df[1, ], ped), "father")
  expect_true(is.na(variant_origin(df[1, ], NULL)))
  df$P002_GT[1] <- "0/0"
  expect_equal(variant_origin(df[1, ], ped), "not in the parents")
  expect_equal(nrow(gene_other_variants(df, df[4, ], ped)), 0)
})
