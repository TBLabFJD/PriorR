# Synthetic demo data for PriorR (user guide, trying PriorR): an invented trio (P001 affected son,
# P002 father, P003 mother) and invented CNVs. Positions, genotypes and annotations are simulated
# (real gene symbols, random positions): they are NOT variants of any person.
# Run from the app/ folder:  Rscript Program_files/scripts/make_demo_data.R
suppressMessages(library(shiny))
for (f in list.files("R", full.names = TRUE)) source(f)
set.seed(2026)
out_dir <- file.path("Program_files", "examples")

# Genes: retinal dystrophy genes (RetiGene) and other genes with constraint data
catalog <- panel_catalog(panels_dir)
retina <- unique(unlist(catalog$genes[catalog$id == "retigene:all" | grepl("^retigene:", catalog$id)]))
others <- setdiff(constraint_table()$SYMBOL, retina)
genes <- c(sample(retina, min(150, length(retina))), sample(others, 650))
chroms <- paste0("chr", c(1:22, "X"))
gene_chrom <- stats::setNames(sample(chroms, length(genes), replace = TRUE, prob = c(rep(1, 22), 0.6)), genes)
gene_start <- stats::setNames(round(runif(length(genes), 1e6, 1.5e8)), genes)
recessive <- stats::setNames(runif(length(genes)) < 0.45, genes)

n <- 4000
gene <- sample(genes, n, replace = TRUE, prob = ifelse(genes %in% retina, 3, 1))
consequences <- c("missense_variant", "synonymous_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
                  "stop_gained", "frameshift_variant", "splice_donor_variant", "splice_region_variant,intron_variant",
                  "inframe_deletion", "upstream_gene_variant", "downstream_gene_variant")
consequence <- sample(consequences, n, replace = TRUE, prob = c(30, 18, 25, 6, 3, 2, 2, 1, 4, 1, 4, 4))
region <- ifelse(grepl("missense|synonymous|stop|frameshift|inframe", consequence), "EXONIC",
          ifelse(grepl("splice_donor|splice_region", consequence), "SPLICING",
          ifelse(grepl("UTR", consequence), ifelse(grepl("3_prime", consequence), "3UTR", "5UTR"),
          ifelse(grepl("upstream", consequence), "UPSTREAM", ifelse(grepl("downstream", consequence), "DOWNSTREAM", "INTRONIC")))))
bases <- c("A", "C", "G", "T")
ref <- sample(bases, n, replace = TRUE)
alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1), "")
indel <- grepl("frameshift|inframe", consequence) | runif(n) < 0.05
ref[indel] <- vapply(which(indel), function(i) paste(sample(bases, sample(2:4, 1), TRUE), collapse = ""), "")
alt[indel] <- substr(ref[indel], 1, 1)
pos <- gene_start[gene] + round(runif(n, 0, 80000))

# Frequencies: most common, some rare, some absent
af <- ifelse(runif(n) < 0.55, runif(n, 0.01, 0.6), ifelse(runif(n) < 0.6, 10^runif(n, -5, -2), NA))
damaging <- grepl("stop|frameshift|splice_donor", consequence)
af[damaging & runif(n) < 0.7] <- NA
num <- function(x, digits = 6) ifelse(is.na(x), NA, signif(x, digits))
cadd <- ifelse(damaging, runif(n, 30, 45), ifelse(consequence == "missense_variant", runif(n, 5, 34), runif(n, 0, 15)))
revel <- ifelse(consequence == "missense_variant", rbeta(n, 1.2, 2.5), NA)
spliceai <- ifelse(grepl("splice", consequence), runif(n, 0, 0.9), runif(n, 0, 0.08))
npred <- ifelse(consequence == "missense_variant", 21, NA)
npath <- ifelse(is.na(revel), NA, round(revel * 21))
clnsig <- ifelse(runif(n) < 0.04, sample(c("Pathogenic", "Likely_pathogenic", "Uncertain_significance", "Benign", "Likely_benign",
                                           "Conflicting_classifications_of_pathogenicity"), n, TRUE, prob = c(2, 2, 4, 3, 3, 1)), NA)
clnsig[damaging & runif(n) < 0.3] <- "Pathogenic"
clnsig[!is.na(af) & af > 0.05 & clnsig %in% c("Pathogenic", "Likely_pathogenic")] <- "Benign"

# Trio genotypes: the son inherits from his parents; a few de novo variants
gt_child <- sample(c("0/1", "1/1"), n, replace = TRUE, prob = c(0.7, 0.3))
from_father <- runif(n) < 0.5
gt_father <- ifelse(gt_child == "1/1", "0/1", ifelse(from_father, sample(c("0/1", "1/1"), n, TRUE, c(0.8, 0.2)), "0/0"))
gt_mother <- ifelse(gt_child == "1/1", "0/1", ifelse(from_father, "0/0", sample(c("0/1", "1/1"), n, TRUE, c(0.8, 0.2))))
# chrX of a son: hemizygous, always from his mother
on_x <- gene_chrom[gene] == "chrX"
gt_child[on_x] <- "1/1"
gt_mother[on_x] <- ifelse(runif(sum(on_x)) < 0.8, "0/1", "1/1")
gt_father[on_x] <- sample(c("0/0", "1/1"), sum(on_x), TRUE)
denovo <- sample(which(gt_child == "0/1" & !on_x & (is.na(af) | af < 1e-4) & region %in% c("EXONIC", "SPLICING")), 6)
gt_father[denovo] <- "0/0"; gt_mother[denovo] <- "0/0"
# a compound heterozygous pair in a recessive retinal gene
ch_gene <- sample(intersect(retina[gene_chrom[retina] != "chrX"], names(which(recessive))), 1)
ch <- sample(which(gene == ch_gene), 2)
gene[ch] <- ch_gene; gt_child[ch] <- "0/1"; gt_father[ch] <- c("0/1", "0/0"); gt_mother[ch] <- c("0/0", "0/1")
consequence[ch] <- c("missense_variant", "stop_gained"); region[ch] <- "EXONIC"; af[ch] <- c(3e-5, NA)
cadd[ch] <- c(29.5, 41); revel[ch] <- c(0.91, NA); clnsig[ch] <- c("Uncertain_significance", "Pathogenic")
depth <- function(k) round(rnorm(n, 35, 9)) |> pmax(8)
genotype_cols <- function(id, gt) {
  dp <- depth()
  alt_reads <- ifelse(gt == "0/0", rbinom(n, dp, 0.003), ifelse(gt == "1/1", dp - rbinom(n, dp, 0.02), rbinom(n, dp, 0.5)))
  stats::setNames(data.frame(gt, round(alt_reads / dp, 2), paste0(dp - alt_reads, "_", alt_reads), dp, stringsAsFactors = FALSE),
                  paste0(id, c("_GT", "_VAF", "_AD", "_DP")))
}
aa <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
cpos <- sample(50:4000, n, TRUE)
hgvsc <- sprintf("NM_%06d.%d:c.%d%s>%s", match(gene, genes) + 1000, sample(1:5, n, TRUE), cpos, ref, alt)
hgvsp <- ifelse(region == "EXONIC", sprintf("NP_%06d.1:p.%s%d%s", match(gene, genes) + 1000, sample(aa, n, TRUE), cpos %/% 3,
                                            ifelse(grepl("stop", consequence), "Ter", ifelse(grepl("synonymous", consequence), "=", sample(aa, n, TRUE)))), NA)
omim <- ifelse(runif(length(genes)) < 0.6, paste0("Simulated disease ", seq_along(genes), ", ",
                                                  ifelse(recessive, "Autosomal recessive", "Autosomal dominant")), NA)
names(omim) <- genes
glow <- stats::setNames(ifelse(genes %in% retina, 0, sample(1:1800, length(genes), TRUE)), genes)

snv <- data.frame(CHROM = gene_chrom[gene], POS = pos, REF = ref, ALT = alt,
                  Location = paste0(gene_chrom[gene], ":", pos), SYMBOL = gene, Gene_full_name = paste(gene, "(synthetic demo)"),
                  GLOWgenes = glow[gene], VARIANT_CLASS = ifelse(indel, ifelse(nchar(ref) > nchar(alt), "deletion", "insertion"), "SNV"),
                  genotype_cols("P001", gt_child), genotype_cols("P002", gt_father), genotype_cols("P003", gt_mother),
                  Existing_variation = NA, Genomic_region = region, CANONICAL = "YES", Feature = sub(":.*", "", hgvsc),
                  BIOTYPE = "protein_coding", Consequence = consequence, HGVSc = hgvsc, HGVSp = hgvsp,
                  CLNSIG = clnsig, OMIM_phenotype = omim[gene], Orphanet_disorder = NA,
                  gnomADg_AF = num(af), gnomADg_AF_popmax = num(ifelse(is.na(af), NA, pmin(1, af * runif(n, 1, 1.8)))),
                  gnomADg_popmax = ifelse(is.na(af), NA, sample(c("nfe", "afr", "amr", "eas", "sas"), n, TRUE)),
                  gnomADg_nhomalt = ifelse(is.na(af), NA, round(af^2 * 76000)),
                  gnomADe_AF = num(ifelse(region %in% c("EXONIC", "SPLICING") & !is.na(af), af * runif(n, 0.8, 1.2), NA)),
                  gnomADe_AF_popmax = num(ifelse(region %in% c("EXONIC", "SPLICING") & !is.na(af), af * runif(n, 1, 1.8), NA)),
                  FJD_MAF_AF = num(ifelse(is.na(af), NA, pmin(1, af * runif(n, 0.5, 1.5)))), CSVS_AF = num(ifelse(is.na(af), NA, af)),
                  CADD_PHRED = round(cadd, 2), MutScore = round(ifelse(is.na(revel), NA, pmin(1, revel + rnorm(n, 0, 0.1))), 3),
                  REVELScore = round(revel, 3), N_Pathogenic_pred = npath, N_Benign_pred = npred - npath, N_predictions = npred,
                  Pathogenic_pred = ifelse(is.na(npath) | npath == 0, NA, "SIFT,PolyPhen,MetaRNN_pred,ClinPred_pred"),
                  Benign_pred = ifelse(is.na(npath) | npath == 21, NA, "LRT_pred,PROVEAN_pred"),
                  SpliceAI_DS_Max = round(spliceai, 2), ada_score = ifelse(grepl("splice", consequence), round(runif(n), 3), NA),
                  rf_score = ifelse(grepl("splice", consequence), round(runif(n), 3), NA), ExACpLI = NA,
                  stringsAsFactors = FALSE, check.names = FALSE)
snv <- snv[order(match(snv$CHROM, chroms), snv$POS), ]
snv <- snv[!duplicated(snv[, c("CHROM", "POS", "REF", "ALT")]), ]
data.table::fwrite(snv, file.path(out_dir, "demo_trio_synthetic.tsv"), sep = "\t", na = "")

# CNVs (AnnotSV-like columns) of P001
m <- 80
cnv_gene <- replicate(m, paste(sample(genes, sample(1:4, 1)), collapse = ";"))
start <- round(runif(m, 1e6, 1.4e8)); len <- round(10^runif(m, 3, 6.5))
type <- sample(c("DEL", "DUP"), m, TRUE)
cnv <- data.frame(AnnotSV_ID = sprintf("SIM_%03d", 1:m), SV_chrom = sub("chr", "", sample(chroms, m, TRUE)), SV_start = start, SV_end = start + len,
                  SV_length = ifelse(type == "DEL", -len, len), SV_type = type, Samples_ID = "P001",
                  P001 = sample(c("0/1", "1/1"), m, TRUE, c(0.85, 0.15)), Annotation_mode = "full",
                  CytoBand = paste0(sample(c("p", "q"), m, TRUE), sample(11:36, m, TRUE)), Gene_name = cnv_gene,
                  Gene_count = lengths(strsplit(cnv_gene, ";")), Location = sample(c("txStart-txEnd", "exon3-exon7", "intron2-exon5"), m, TRUE),
                  Location2 = sample(c("CDS", "UTR", "CDS;UTR"), m, TRUE),
                  B_gain_AFmax = num(ifelse(type == "DUP" & runif(m) < 0.4, runif(m, 0, 0.1), NA)),
                  B_loss_AFmax = num(ifelse(type == "DEL" & runif(m) < 0.4, runif(m, 0, 0.1), NA)),
                  OMIM_phenotype = omim[sub(";.*", "", cnv_gene)], AnnotSV_ranking_score = round(runif(m, -1, 1.2), 2),
                  ACMG_class = sample(1:5, m, TRUE, c(3, 3, 5, 2, 1)), stringsAsFactors = FALSE)
data.table::fwrite(cnv, file.path(out_dir, "demo_cnv_synthetic.tsv"), sep = "\t", na = "")
message("Demo data: ", nrow(snv), " SNVs (", length(denovo), " de novo, compound heterozygous in ", ch_gene, "), ", m, " CNVs")
