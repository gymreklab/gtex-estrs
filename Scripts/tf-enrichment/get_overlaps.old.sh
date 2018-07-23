#!/bin/bash

TISSUE=$1
INDEXDIR=$2
PREFIX=$3

source params.sh
 
# Get eSTR/all STRs    Causal.eSTRs / all STRs
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Master.table | \
#	awk '($11<0.1)' | awk '($11>0)' |awk '$7 ~ /STR_/ ' | awk '{print $2 "\t" $3 "\t" $3+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_Causal.bed
	 awk '($11<0.1)'| awk '($11>0)' | awk '{print $2 "\t" $3 "\t" $3+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed
#    awk '($4<0.1)'| awk '$10 ~ /STR_/ ' | awk '{print $2 "\t" $4 "\t" $4+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_Causal_eSTRs.bed #chr - str.start - str.start+1


#   
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Lin_Reg_Out | \
	grep -v gene | awk '{print $2 "\t" $4 "\t" $4+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed	
#    grep -v gene | awk '{print $3 "\t" $5 "\t" $5+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed #chrom - str.start - str.start+1


#bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_Causal.bed
bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed
bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed



# Giggle search each list  [Change eSTRs to Causal_eSTRs
giggle search -q ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_overlap.tab      #_Causal_overlap.tab
giggle search -q ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_${PREFIX}_allSTR_overlap.tab


# Join tables to get numbers for Fisher exact test
# (num_estr_tf, num_notestr_tf, num_estr_nottf, num_notestr_nottf) 
#eSTRs  to Causal_eSTRs
totalestr=$(zcat ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed | wc -l)  #_Causal.bed | wc -l)
totalstr=$(zcat ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed | wc -l)
join -t $'\t' ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_overlap.tab ${TMPDIR}/${TISSUE}_${PREFIX}_allSTR_overlap.tab | \
    awk -v"totalestr=$totalestr" -v "totalstr=$totalstr" \
    '{print $1 "\t" $2 "\t" $3-$2 "\t" (totalestr-$2) "\t" (totalstr-($3-$2+totalestr))}' > ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_table.tab  #${TMPDIR}/${TISSUE}_${PREFIX}_Causal_table.tab



# Perform Fisher exact test
./fisher_exact.py ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_table.tab > ${OUTDIR}/${TISSUE}_${PREFIX}_eSTRs_enrich.tab  #${OUTDIR}/${TISSUE}_${PREFIX}_Causal_enrich.tab

#TODO
# REPEAT on CAUSAL eSTRS  to output ${OUTDIR}/${TISSUE}_${PREFIX}_enrich_causal.tab - TODO 
# join ${OUTDIR}/${TISSUE}_${PREFIX}_enrich_causal.tab  ${OUTDIR}/${TISSUE}_${PREFIX}_enrich.tab, output category and odds_ratio_causal/odds_ratio_allestr - TODO
