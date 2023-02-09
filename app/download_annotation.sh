
# dockershareddir="/home/gonzalo/tblab/mnt/tblab/gonzalo/pruebas_vep_priorr/priorr_shared_folder"
# docker run -it -v ${dockershareddir}:${dockershareddir} -u $(id -u):$(id -g) priorr_annot:latest bash

# ######################
# # Para el dockerfile #
# ######################

# apt-get update && apt-get -y install \
# build-essential \
# cpanminus \
# curl \
# libmysqlclient-dev \
# libpng-dev \
# libssl-dev \
# zlib1g-dev \
# libbz2-dev \
# liblzma-dev \
# locales \
# openssl \
# perl \
# perl-base \
# unzip \
# vim \
# git && \
# apt-get -y purge manpages-dev && \
# apt-get clean && \
# rm -rf /var/lib/apt/lists/* # buildkit






# bcftools
# bgzip
# tabix
# vep
# automap
# R
# R (optparse)








######################
# VEP cache download #
######################

dockershareddir=$1
assembly=$2 # GRCh37 or GRCh38 or both

# dockershareddir="/mnt/tblab/gonzalo/pruebas_vep_priorr/priorr_shared_folder"
# assembly="GRCh37" # GRCh37 or GRCh38 or both



if [[ "$assembly" == "GRCh37" || "$assembly" == "both" ]]; then 
perl INSTALL.pl --AUTO cfp \
--NO_UPDATE \
--SPECIES homo_sapiens_refseq \
--ASSEMBLY GRCh37 \
--CACHEDIR ${dockershareddir}/.vep \
--CACHE_VERSION 105 \
-g dbscSNV,LoFtool,ExACpLI,dbNSFP,MaxEntScan,CADD
zcat ${dockershareddir}/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz > ${dockershareddir}/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
bgzip --force ${dockershareddir}/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
fi


if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then 
perl INSTALL.pl --AUTO cfp \
--NO_UPDATE \
--SPECIES homo_sapiens_refseq \
--ASSEMBLY GRCh38 \
--CACHEDIR ${dockershareddir}/.vep \
--CACHE_VERSION 105 \
-g dbscSNV,LoFtool,ExACpLI,dbNSFP,MaxEntScan,CADD
zcat ${dockershareddir}/.vep/homo_sapiens_refseq/105_GRCh28/Homo_sapiens.GRCh28.75.dna.primary_assembly.fa.gz > ${dockershareddir}/.vep/homo_sapiens_refseq/105_GRCh28/Homo_sapiens.GRCh28.75.dna.primary_assembly.fa
bgzip --force ${dockershareddir}/.vep/homo_sapiens_refseq/105_GRCh28/Homo_sapiens.GRCh28.75.dna.primary_assembly.fa

fi



mkdir ${dockershareddir}/custom
cd ${dockershareddir}/custom

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
# rm -r ./tmp_sort/*
fi

if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then
mkdir ./tmp_sort/ 
zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | sort -T ./tmp_sort/ -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.3a_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.3a_grch38.gz
# rm -r ./tmp_sort/*
fi

zcat dbNSFP4.3_gene.complete.gz | cut -f1,14,20,26,27,30,31,32,33,36,37,95 > dbNSFP4.3_gene.complete.pvm.txt

# rm dbNSFP4.3a_variant.chr* tryhg* try.vcf search_dbNSFP43a* LICENSE.txt h dbNSFP4.3a.zip dbNSFP4.3_gene.gz dbNSFP4.3_gene.complete.gz



# dbscSNV (April 12, 2015)
# dbscSNV includes all potential human SNVs within splicing consensus regions (−3 to +8 at the 5’ splice site and −12 to +2 at the 3’ splice site), i.e. scSNVs, related functional annotations and two ensemble prediction scores for predicting their potential of altering splicing.
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip
unzip dbscSNV1.1.zip
head -n1 dbscSNV1.1.chr1 > h

if [[ "$assembly" == "GRCh37" || "$assembly" == "both" ]]; then 
cat dbscSNV1.1.chr* | grep -v ^chr | cat h - | bgzip -c > dbscSNV1.1_GRCh37.txt.gz
tabix -s 1 -b 2 -e 2 -c c dbscSNV1.1_GRCh37.txt.gz
fi

if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then 
cat dbscSNV1.1.chr* | grep -v ^chr | sort -k5,5 -k6,6n | cat h - | awk '$5 != "."' | bgzip -c > dbscSNV1.1_GRCh38.txt.gz
tabix -s 5 -b 6 -e 6 -c c dbscSNV1.1_GRCh38.txt.gz
fi

rm dbscSNV1.1.chr* h dbscSNV1.1.zip




# gnomAD
if [[ "$assembly" == "GRCh37" || "$assembly" == "both" ]]; then 
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi

bcftools annotate -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,\
INFO/popmax,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/nhomalt_popmax \
--threads 8 -O z -o gnomad.exomes.r2.1.1.sites.annotfilt.vcf.bgz gnomad.exomes.r2.1.1.sites.vcf.bgz
tabix -p vcf gnomad.exomes.r2.1.1.sites.annotfilt.vcf.bgz

rm gnomad.exomes.r2.1.1.sites.vcf.bgz gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
fi


if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then 
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi

bcftools annotate -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,\
INFO/popmax,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/nhomalt_popmax \
--threads 8 -O z -o gnomad.exomes.r2.1.1.sites.liftover_grch38.annotfilt.vcf.bgz gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
tabix -p vcf gnomad.exomes.r2.1.1.sites.liftover_grch38.annotfilt.vcf.bgz

rm gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi

fi



# ExACpLI
# A VEP plugin that adds the probabililty of a gene being loss-of-function intolerant (pLI) to the VEP output.
# The closer pLI is to 1, the more likely the gene is loss-of-function (LoF) intolerant.
# wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/ExACpLI_values.txt
# mv ExACpLI_values.txt ${VEP_CACHE}/Plugins/


# LoFtool
# LoFtool provides a rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations for each gene
# The lower the LoFtool gene score percentile the most intolerant is the gene to functional variation.
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/LoFtool_scores.txt
# mv LoFtool_scores.txt ${VEP_CACHE}/Plugins/


# MaxEntScan
# 
# http://genes.mit.edu/burgelab/maxent/download/
wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
gunzip fordownload.tar.gz
tar -xvf fordownload.tar
mv fordownload maxEntScan
rm fordownload.tar



# ClinVar
if [[ "$assembly" == "GRCh37" || "$assembly" == "both" ]]; then 
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
mv clinvar.vcf.gz clinvar.GRCh37.vcf.gz
mv clinvar.vcf.gz.tbi clinvar.GRCh37.vcf.gz.tbi
fi

if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then 
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
mv clinvar.vcf.gz clinvar.GRCh38.vcf.gz
mv clinvar.vcf.gz.tbi clinvar.GRCh38.vcf.gz.tbi
fi



# MutScore
if [[ "$assembly" == "GRCh37" || "$assembly" == "both" ]]; then 
wget https://storage.googleapis.com/rivolta_mutscore/mutscore-v1.0-hg19.tsv.gz
vcfFromBed="mutscore-v1.0-hg19.vcf.gz"
bcftools view -h clinvar.GRCh37.vcf.gz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h clinvar.GRCh37.vcf.gz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=Score,Number=A,Type=String,Description="MutScore">' | bgzip -c  >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c >> ${vcfFromBed}
zcat mutscore-v1.0-hg19.tsv.gz | awk -F"\t" '{print "chr"$1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tScore="$5}' | bgzip -c >> ${vcfFromBed}
tabix -p vcf ${vcfFromBed}
rm mutscore-v1.0-hg19.tsv.gz*
fi

if [[ "$assembly" == "GRCh38" || "$assembly" == "both" ]]; then 
wget https://storage.googleapis.com/rivolta_mutscore/mutscore-v1.0-hg38.tsv.gz
vcfFromBed="mutscore-v1.0-hg38.vcf.gz"
bcftools view -h clinvar.GRCh38.vcf.gz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h clinvar.GRCh38.vcf.gz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=Score,Number=A,Type=String,Description="MutScore">' | bgzip -c  >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c >> ${vcfFromBed}
zcat mutscore-v1.0-hg19.tsv.gz | awk -F"\t" '{print "chr"$1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tScore="$5}' | bgzip -c >> ${vcfFromBed}
tabix -p vcf ${vcfFromBed}
rm mutscore-v1.0-hg38.tsv.gz*
fi



