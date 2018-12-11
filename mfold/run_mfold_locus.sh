#!/bin/bash

# Perform mfold analysis on one locus
# Usage:
# ./run_mfold_locus.sh <chrom:start-end>
# 
# Generates file <chrom_start>.energies.tab
# Uses SSC VCFs to get alleles. Not ideal, but ALT field of gtex VCF is weird
# Should probably modify to use GTEx VCFs

source params.sh

# Inputs
LOCUS=$1 # 21:45196326-45196326
LOCNAME=$(echo $LOCUS | cut -f 1-2 -d'_' | sed 's/:/_/') 

# Get VCF file to use
chrom=$(echo $LOCUS | cut -f 1 -d':')
VCF=${VCFDIR}/hipstr.chr${chrom}.allfilters.vcf.gz

# Make temporary directory
TMPDIR=$(mktemp --tmpdir=/tmp -d mfoldXXXXXXX)
echo $TMPDIR
# Extract flanking regions and alleles for this locus
seq=$(bcftools query -r ${LOCUS} -f "%CHROM\t%POS\t%REF\n" $VCF | \
    awk -v"window=$WINDOW" '{print $1 "\t" $2-window-1 "\t" $2+length($3)+window-1}' | \
    bedtools getfasta -bed stdin -fi ${REFFA} | tail -n 1)
lflank=$(echo $seq | awk -v"window=$WINDOW" '{print substr($1, 1, window)}')
rflank=$(echo $seq | awk -v"window=$WINDOW" '{print substr($1, length($1)-window+1, window)}')

echo $LOCUS

# Compute free energy for each allele
cd ${TMPDIR}
alleles=$(extract_common_alleles.py $VCF $LOCUS $MINCOUNT)
#alleles=$(bcftools query -r ${LOCUS} -f "%REF,%ALT\n" $VCF | sed 's/,/\t/g')
TMPFA=test_seqs.fa
for a in $alleles
do
    # Forward direction
#    echo $a | awk '{print ">" length($0) ":" $0}' > ${TMPFA}
#    echo ${lflank}${a}${rflank} >> ${TMPFA}
#    mfold_MG SEQ=${TMPFA} NA=DNA > /dev/null 2>&1
#    head -n 1 test_seqs.fa.ct > ${TMPDIR}/${a}_dnaplus.ct
#    rm ${TMPFA}*
    echo ${a} | dna2rna |  awk '{print ">" length($0) ":" $0}' > ${TMPFA}
    echo ${lflank}${a}${rflank} | dna2rna >> ${TMPFA}
    mfold_MG SEQ=${TMPFA} NA=RNA > /dev/null 2>&1
    head -n 1 test_seqs.fa.ct > ${TMPDIR}/${a}_rnaplus.ct
    rm ${TMPFA}*
    # Reverse direction
#    echo $a | reverse_complement | awk '{print ">" length($0) ":" $0}' > ${TMPFA}
#    echo ${lflank}${a}${rflank} | reverse_complement >> ${TMPFA}
#    mfold_MG SEQ=${TMPFA} NA=DNA > /dev/null 2>&1
#    head -n 1 test_seqs.fa.ct > ${TMPDIR}/${a}_dnaminus.ct
#    rm ${TMPFA}*
    echo ${a} | reverse_complement | dna2rna | awk '{print ">" length($0) ":" $0}' > ${TMPFA}
    echo ${lflank}${a}${rflank} | reverse_complement | dna2rna>> ${TMPFA}    
    mfold_MG SEQ=${TMPFA} NA=RNA > /dev/null 2>&1
    head -n 1 test_seqs.fa.ct > ${TMPDIR}/${a}_rnaminus.ct
    rm ${TMPFA}*
done

#OUTFILE=${OUTDIR}/perlocus/${chrom}/${LOCNAME}.energies.tab
#cat ${TMPDIR}/*_dnaplus.ct | awk -v"name=$LOCNAME" '{split($NF,a,":"); print name "\tDNA\t+\t" a[1]"\t" $NF "\t" $4}' > ${OUTFILE}
#cat ${TMPDIR}/*_dnaminus.ct | awk -v"name=$LOCNAME" '{split($NF,a,":"); print name "\tDNA\t-\t" a[1]"\t" $NF "\t" $4}' >> ${OUTFILE}

OUTFILE=${OUTDIR}/perlocus/${chrom}/${LOCNAME}.rna.energies.tab
cat ${TMPDIR}/*_rnaplus.ct | awk -v"name=$LOCNAME" '{split($NF,a,":"); print name "\tRNA\t+\t" a[1]"\t" $NF "\t" $4}' > ${OUTFILE}
cat ${TMPDIR}/*_rnaminus.ct | awk -v"name=$LOCNAME" '{split($NF,a,":"); print name "\tRNA\t-\t" a[1]"\t" $NF "\t" $4}' >> ${OUTFILE}


#rm -rf ${TMPDIR}
