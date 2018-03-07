#!/bin/bash

TISSUE=$1

source params.sh

# Get eSTR/all STRs
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/PQValues | \
    awk '($8<0.1)' | awk '{print $2 "\t" $4 "\t" $4+1}' > ${TMPDIR}/${TISSUE}_eSTRs.bed
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Lin_Reg_Out | \
    grep -v gene | awk '{print $3 "\t" $5 "\t" $5+1}' > ${TMPDIR}/${TISSUE}_allSTRs.bed
bgzip -f ${TMPDIR}/${TISSUE}_eSTRs.bed
bgzip -f ${TMPDIR}/${TISSUE}_allSTRs.bed

# Giggle search each list
giggle search -q ${TMPDIR}/${TISSUE}_eSTRs.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_eSTR_overlap.tab
giggle search -q ${TMPDIR}/${TISSUE}_allSTRs.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_allSTR_overlap.tab


