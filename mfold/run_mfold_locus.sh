#!/bin/bash

# Perform mfold analysis on one locus
# Usage:
# ./run_mfold_locus.sh <chrom:start-end>
# 
# Generates file <chrom_start>.energies.tab
# Uses SSC VCFs to get alleles. Not ideal, but ALT field of gtex VCF is weird
# Should probably modify to use GTEx VCFs

# Inputs
LOCUS=$1 # 21:45196326-45196326
LOCNAME=$(echo $LOCUS | cut -f 1-2 -d'_' | sed 's/:/_/') 

# Other parameters 
WINDOW=50 # Use this much context around the STR
REFFA=/storage/resources/dbase/human/hs37d5/hs37d5.fa

# Get VCF file to use
chrom=$(echo $LOCUS | cut -f 1 -d':')
VCF=/storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr${chrom}.allfilters.vcf.gz

# Make temporary directory
TMPDIR=$(mktemp -d)

# Extract flanking regions and alleles for this locus
seq=$(bcftools query -r ${LOCUS} -f "%CHROM\t%POS\t%REF\n" $VCF | \
    awk -v"window=$WINDOW" '{print $1 "\t" $2-window-1 "\t" $2+length($3)+window-1}' | \
    bedtools getfasta -bed stdin -fi ${REFFA} | tail -n 1)
lflank=$(echo $seq | awk -v"window=$WINDOW" '{print substr($1, 1, window)}')
rflank=$(echo $seq | awk -v"window=$WINDOW" '{print substr($1, length($1)-window+1, window)}')

# Compute free energy for each allele
cd ${TMPDIR}
alleles=$(bcftools query -r ${LOCUS} -f "%REF,%ALT\n" $VCF | sed 's/,/\t/g')
for a in $alleles
do
    TMPFA=test_seqs.fa
    echo ">"${a} > ${TMPFA}
    echo ${lflank}${a}${rflank} >> ${TMPFA}
    mfold SEQ=${TMPFA} NA=DNA 
    head -n 1 test_seqs.fa.ct > ${TMPDIR}/${a}.ct
    rm test_seqs.fa*
done
cd -

cat ${TMPDIR}/*.ct | awk -v"name=$LOCNAME" '{print name "\t" length($NF) "\t" $NF "\t" $4}' > ${LOCNAME}.energies.tab
