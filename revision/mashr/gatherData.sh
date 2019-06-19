#!/bin/bash

#example where I learned this from
# ls -ld /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/* | grep '^d' | sed -e 's/\s\+/ /g' | cut -f 9 -d ' '

tissuetypes='Adipose-Subcutaneous
Adipose-Visceral
Artery-Aorta
Artery-Tibial
Brain-Caudate
Brain-Cerebellum
Cells-Transformedfibroblasts
Esophagus-Mucosa
Esophagus-Muscularis
Heart-LeftVentricle
Lung
Muscle-Skeletal
Nerve-Tibial
Skin-NotSunExposed
Skin-SunExposed
Thyroid
WholeBlood
'
#Don't seem to have tables
#'Kidney-Cortex Liver '

#escape the parentheses in the tissue types names
tissuetypes=$(printf '%q ' $tissuetypes)

workdir='/storage/mgymrek/gtex-estrs/revision/mashr/'
mkdir -p ${workdir}/output-strs
mkdir -p ${workdir}/input-strs
mkdir -p ${workdir}/intermediate-strs
mkdir -p ${workdir}/output-snps
mkdir -p ${workdir}/input-snps
mkdir -p ${workdir}/intermediate-snps

datadir=/storage/mgymrek/gtex-estrs/revision/
### Get STR data ###
for tissue in $tissuetypes; do
    echo "gene,chrom,str.start,beta,significant,beta.se" | sed 's/,/\t/g' > ${workdir}/input-strs/$tissue.table
    cat ${datadir}/strreg/${tissue}_strreg.tab | grep -v chrom | \
	awk -F"\t" '{print $2 "\t" $3 "\t" $5 "\t" $10 "\t" ($13<10**-4) "\t" $11}' >> ${workdir}/input-strs/$tissue.table
done

### Get SNP data ###
for tissue in $tissuetypes; do
    echo "gene,chrom,str.start,beta,significant,beta.se" | sed 's/,/\t/g' > ${workdir}/input-snps/$tissue.table
    cat ${datadir}/snpreg/${tissue}_snpreg.tab | grep -v chrom | \
	awk -F"\t" '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" ($9<10**-4) "\t" $7}' >> ${workdir}/input-snps/$tissue.table
done
