#!/bin/bash


dockershareddir=$1
assembly=$2 # GRCh37 or GRCh38 or both

mkdir ${dockershareddir}/custom
cd ${dockershareddir}/custom


#share_dir='/home/raquel/Priorr_test/custom'
#cd $share_dir

# dbNSFP (February 18, 2022)
# dbNSFP is a database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome.
# It compiles prediction scores from 38 prediction algorithms (SIFT, SIFT4G, Polyphen2-HDIV, Polyphen2-HVAR, LRT, MutationTaster2, MutationAssessor, FATHMM, MetaSVM, MetaLR, MetaRNN, CADD, CADD_hg19, VEST4, PROVEAN, FATHMM-MKL coding, FATHMM-XF coding, fitCons x 4, LINSIGHT, DANN, GenoCanyon, Eigen, Eigen-PC, M-CAP, REVEL, MutPred, MVP, MPC, PrimateAI, GEOGEN2, BayesDel_addAF, BayesDel_noAF, ClinPred, LIST-S2, ALoFT), 9 conservation scores (PhyloP x 3, phastCons x 3, GERP++, SiPhy and bStatistic) and other related information including allele frequencies observed in the 1000 Genomes Project phase 3 data, UK10K cohorts data, ExAC consortium data, gnomAD data and the NHLBI Exome Sequencing Project ESP6500 data, various gene IDs from different databases, functional descriptions of genes, gene expression and gene interaction information, etc.
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.3a.zip
unzip dbNSFP4.3a.zip
zcat dbNSFP4.3a_variant.chr1.gz | head -n1 > h



if [[ "$assembly" == "GRCh37" || "$assembly" == "both" ]]; then
mkdir ./tmp_sort/
zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | awk '$8 != "." ' | sort -T ./tmp_sort/ -k8,8 -k9,9n - | cat h - | bgzip -c > dbNSFP4.3a_grch37.gz
tabix -s 8 -b 9 -e 9 dbNSFP4.3a_grch37.gz
#  rm -r ./tmp_sort/*
fi

if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then
mkdir ./tmp_sort/
zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | sort -T ./tmp_sort/ -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.3a_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.3a_grch38.gz
# rm -r ./tmp_sort/*
fi

zcat dbNSFP4.3_gene.complete.gz | cut -f1,14,20,26,27,30,31,32,33,36,37,95 > dbNSFP4.3_gene.complete.pvm.txt

