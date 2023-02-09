#!/bin/bash



#docker run -it -v /mnt/tblab/gonzalo/pruebas_vep_priorr/:/mnt/tblab/gonzalo/pruebas_vep_priorr/ \
#-u $(id -u):$(id -g) priorr_gonzalo bash


sharedir=$1
assembly=$2 # GRCh37 or GRCh38
inputvcf=$3

#cd /mnt/tblab/gonzalo/pruebas_vep_priorr/example_dir/


#inputvcf="/mnt/tblab/gonzalo/pruebas_vep_priorr/example_dir/22-2341.final.gatk.vcf"
#sharedir="/mnt/tblab/gonzalo/pruebas_vep_priorr/priorr_shared_folder/"


vep_cache="${sharedir}/.vep"
vep_plugins="${sharedir}/.vep/Plugins"
loFtool="${sharedir}/custom/LoFtool_scores.txt"
exACpLI="${sharedir}/custom/ExACpLI_values.txt"
maxEntScan="${sharedir}/custom/maxEntScan/"
dbNSFP_gene="${sharedir}/custom/dbNSFP4.3_gene.complete.pvm.txt"

if [ "$assembly" == 'GRCh37' ] ; then
automap_assembly="hg19"
vep_fasta="${sharedir}/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
dbNSFP="${sharedir}/custom/dbNSFP4.3a_grch37.gz"
dbscSNV="${sharedir}/custom/dbscSNV1.1_GRCh37.txt.gz"
mutScore="${sharedir}/custom/mutscore-v1.0-hg19.vcf.gz"
cLINVAR="${sharedir}/custom/clinvar.GRCh37.vcf.gz"
gNOMADe="gnomad.exomes.r2.1.1.sites.annotfilt.vcf.bgz"

elif [ "$assembly" == 'GRCh38' ] ; then 
automap_assembly="hg38"
vep_fasta="${sharedir}/.vep/homo_sapiens_refseq/105_GRCh38/Homo_sapiens.GRCh38.75.dna.primary_assembly.fa.gz"
dbNSFP="${sharedir}/custom/dbNSFP4.3a_grch38.gz"
dbscSNV="${sharedir}/custom/dbscSNV1.1_GRCh38.txt.gz"
mutScore="${sharedir}/custom/mutscore-v1.0-hg38.vcf.gz"
cLINVAR="${sharedir}/custom/clinvar.GRCh38.vcf.gz"
gNOMADe="gnomad.exomes.r2.1.1.sites.liftover_grch38.annotfilt.vcf.bgz"

fi


maf="1"
regiondict="/home/app/dict_region.csv"




sample="$(basename ${inputvcf} | sed 's/\..*//')"


#######################
# INPUT NORMALIZATION #
#######################

bcftools norm -c s -O z -o ${sample}.normalized.vcf.gz -f ${vep_fasta} ${inputvcf}

inputvcf="${sample}.normalized.vcf.gz"


##########################################
# SAMPLE INFORMATION PREPARATION 4 ANNOT #
##########################################

bcftools view -h ${inputvcf} | grep "##" > ${sample}.vcf_to_annotate.vcf
echo "##INFO=<ID=variant_id,Number=.,Type=String,Description=\"variant identification\">" >> ${sample}.vcf_to_annotate.vcf
echo "##INFO=<ID=Original_pos,Number=.,Type=String,Description=\"original position\">" >> ${sample}.vcf_to_annotate.vcf
for muestra in $(bcftools query -l ${inputvcf})
do
echo "##INFO=<ID=${muestra}_GT,Number=.,Type=String,Description=\"${muestra} Genotype\">" >> ${sample}.vcf_to_annotate.vcf
echo "##INFO=<ID=${muestra}_AD,Number=.,Type=String,Description=\"${muestra} Allelic depths for the ref and alt alleles in the order listed\">" >> ${sample}.vcf_to_annotate.vcf
echo "##INFO=<ID=${muestra}_DP,Number=.,Type=String,Description=\"${muestra} Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">" >> ${sample}.vcf_to_annotate.vcf
echo "##INFO=<ID=${muestra}_GQ,Number=.,Type=String,Description=\"${muestra} Genotype Quality\">" >> ${sample}.vcf_to_annotate.vcf
done
bcftools view -h ${inputvcf} | grep "#CHROM" | cut -f1-8 >> ${sample}.vcf_to_annotate.vcf


bcftools query -f 'variant_id=%CHROM\_%POS\_%REF\_%ALT;Original_pos=%POS;[;%SAMPLE\_GT=%GT][;%SAMPLE\_AD=%AD][;%SAMPLE\_DP=%DP][;%SAMPLE\_GQ=%GQ]\n' ${inputvcf} | sed 's/;//2' | sed 's/,/_/g' > new_info.txt
bcftools view -H ${inputvcf} | cut -f1-8 > old_info.txt
paste -d ';' old_info.txt new_info.txt >> ${sample}.vcf_to_annotate.vcf
rm old_info.txt new_info.txt

bgzip ${sample}.vcf_to_annotate.vcf
tabix -p vcf ${sample}.vcf_to_annotate.vcf.gz


if bcftools view -h ${inputvcf} | grep -Fq hiConfDeNovo; then
fields=",hiConfDeNovo,loConfDeNovo,variant_id,Original_pos"
else
fields=",variant_id,Original_pos"
fi

for muestra in $(bcftools query -l ${inputvcf}); do fields="$(echo "${fields},${muestra}_GT,${muestra}_AD,${muestra}_DP,${muestra}_GQ")"; done
# echo ${fields} > ${sample}.fields.txt



########################
# CANONICAL CHR FILTER #
########################

/ensembl-vep-release-105/filter_vep \
-i ${inputvcf} -o ${sample}.canonicalchr.vcf \
--filter "(CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY)" \
--force_overwrite




#######
# VEP #
#######

/ensembl-vep-release-105/vep \
--cache --offline --dir_cache ${vep_cache} --dir_plugins ${vep_plugins} \
--refseq --species homo_sapiens --assembly ${vep_assembly} --force_overwrite --use_transcript_ref \
--verbose --fork $(nproc) --tab --format vcf --no_stats \
--fasta ${vep_fasta} \
--input_file ${sample}.canonicalchr.vcf \
--output_file ${sample}.vep.tsv \
--check_existing --canonical --numbers --hgvs --biotype --regulatory --symbol --protein \
--sift p --polyphen p --allele_number --variant_class --pubmed \
--plugin dbNSFP,${dbNSFP},\
LRT_pred,M-CAP_pred,MetaLR_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,\
FATHMM_pred,MetaRNN_pred,PrimateAI_pred,DEOGEN2_pred,BayesDel_addAF_pred,BayesDel_noAF_pred,ClinPred_pred,\
LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,\
phyloP30way_mammalian,phastCons30way_mammalian,GERP++_RS,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue,\
CADD_raw,CADD_phred \
--plugin dbscSNV,${dbscSNV} \
--plugin LoFtool,${loFtool} \
--plugin ExACpLI,${exACpLI} \
--plugin MaxEntScan,${maxEntScan} \
--custom ${gNOMADe},gnomADe,vcf,exact,0,AF,AC,AN,nhomalt,popmax,AF_popmax,AC_popmax  \
--custom ${mutScore},Mut,vcf,exact,0,Score \
--custom ${cLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
--custom ${sample}.vcf_to_annotate.vcf.gz,SAMPLE,vcf,exact,0${fields} 






###########
# Automap #
###########

if [[ $(bcftools query -l ${inputvcf} | wc -l) -gt 1 ]]; then
for sample in $(bcftools query -l ${inputvcf}); do

bcftools view -s ${sample} -O v -o ${sample}.indv.vcf ${inputvcf}

if [[ $(bcftools view -H ${sample}.indv.vcf | wc -l) -gt 10000 ]]; then

/home/docker/AutoMap/AutoMap_v1.2.sh \
--vcf ${sample}.indv.vcf \
--out . \
--genome ${automap_assembly}

mv ${sample}/* .
fi
done

else
if [[ $(bcftools view -H ${inputvcf} | wc -l) -gt 10000 ]]; then

/home/docker/AutoMap/AutoMap_v1.2.sh \
--vcf ${inputvcf} \
--out . \
--genome ${automap_assembly}

mv */* .
fi
fi








#########################
# POST-VEP MODIFICATION #
#########################

header_row="$(head -n 1000 ${sample}.vep.tsv | grep "#Uploaded_variation" -n | sed 's/:.*//')"

Rscript /home/app/post-VEP_modification.R \
--input ${sample}.vep.tsv \
--output ${sample}.pvm.tsv \
--numheader ${header_row} \
--dbNSFPgene ${dbNSFP_gene} \
--regiondict ${regiondict} \
--automap ./ \
--maf ${maf}
