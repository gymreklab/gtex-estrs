#!/bin/bash

TISSUE=$1
INDEXDIR=$2
PREFIX=$3

source params.sh

			###############################################
        		####            DATA COLLECT            #######
        		###############################################

###############   CAUSAL eSTRs	################
 #[take chrom  start end/start+1 where best.variariant is STR]
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Master.table | \
	awk '($8<0.1)' | awk '($8>0)' |awk '$6 ~ /STR_/ ' | awk '{print $2 "\t" $3 "\t" $3+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_Causal.bed
#compress
bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_Causal.bed
# Giggle search  list of causal eSTRs
giggle search -q ${TMPDIR}/${TISSUE}_${PREFIX}_Causal.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_${PREFIX}_Causal_overlap.tab

###############   just  eSTRs 	################
#[take chrom  start end/start+1 where qval<0.1 & >0]
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Master.table | \
	awk '($8<0.1)'| awk '($8>0)' | awk '{print $2 "\t" $3 "\t" $3+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed
#compress
bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed
# Giggle search  list of  eSTRs
giggle search -q ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_overlap.tab


#############   Non causal eSTRs ################
#[take chrom  start end/start+1 where qval<0.1 & >0]
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Master.table | \
        awk '($8<0.1)'| awk '($8>0)'|awk '$6 ~ /SNP_/ ' | awk '{print $2 "\t" $3 "\t" $3+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal.bed
#compress
bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal.bed
# Giggle search  list of  eSTRs
giggle search -q ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal_overlap.tab


###############   All    STRs   ################
#[take chrom gene start end/start+1 for all STR]
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Lin_Reg_Out | \
	grep -v gene | awk '{print $2 "\t" $4 "\t" $4+1}' > ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed	
#compress
bgzip -f ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed
# Giggle search  list of STRs
giggle search -q ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed.gz -i ${INDEXDIR} | \
    awk -F "/" '{print $NF}' | sed 's/.bed.gz//' | \
    cut -f 1,3 | sed 's/overlaps://' | sort > ${TMPDIR}/${TISSUE}_${PREFIX}_allSTR_overlap.tab

		###############################################
		####		FISHER TEST		#######
		###############################################

# Making contingency table with numbers for Fisher exact test
#Tables (num_estr_tf, num_notestr_tf, num_estr_nottf, num_notestr_nottf) 
totalestr=$(zcat ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs.bed | wc -l)  
totalstr=$(zcat ${TMPDIR}/${TISSUE}_${PREFIX}_allSTRs.bed | wc -l)
totcausal=$(zcat ${TMPDIR}/${TISSUE}_${PREFIX}_Causal.bed | wc -l)
totncausal=$(zcat ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal.bed | wc -l)

#eSTRs at TFs (All STRs)
join -t $'\t' ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_overlap.tab ${TMPDIR}/${TISSUE}_${PREFIX}_allSTR_overlap.tab | \
    awk -v"totalestr=$totalestr" -v "totalstr=$totalstr" \
    '{print $1 "\t" $2 "\t" $3-$2 "\t" (totalestr-$2) "\t" (totalstr-($3-$2+totalestr))}' > ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_table.tab 
#eCausal at TFs (All STRs)
join -t $'\t' ${TMPDIR}/${TISSUE}_${PREFIX}_Causal_overlap.tab ${TMPDIR}/${TISSUE}_${PREFIX}_allSTR_overlap.tab | \
    awk -v"totcausal=$totcausal" -v "totalstr=$totalstr" \
    '{print $1 "\t" $2 "\t" $3-$2 "\t" (totcausal-$2) "\t" (totalstr-($3-$2+totcausal))}' > ${TMPDIR}/${TISSUE}_${PREFIX}_Causal_table.tab
#NonCausal at TFs (All STRs)
join -t $'\t' ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal_overlap.tab ${TMPDIR}/${TISSUE}_${PREFIX}_allSTR_overlap.tab | \
    awk -v"totncausal=$totncausal" -v "totalstr=$totalstr" \
    '{print $1 "\t" $2 "\t" $3-$2 "\t" (totncausal-$2) "\t" (totalstr-($3-$2+totncausal))}' > ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal_table.tab 

# Perform Fisher exact test
./fisher_exact.py ${TMPDIR}/${TISSUE}_${PREFIX}_eSTRs_table.tab > ${OUTDIR}/${TISSUE}_${PREFIX}_eSTRs_enrich.tab
./fisher_exact.py ${TMPDIR}/${TISSUE}_${PREFIX}_Causal_table.tab > ${OUTDIR}/${TISSUE}_${PREFIX}_Causal_enrich.tab
./fisher_exact.py ${TMPDIR}/${TISSUE}_${PREFIX}_NonCausal_table.tab > ${OUTDIR}/${TISSUE}_${PREFIX}_NonCausal_enrich.tab

#TODO 
# join ${OUTDIR}/${TISSUE}_${PREFIX}_enrich_causal.tab  ${OUTDIR}/${TISSUE}_${PREFIX}_enrich.tab, output category and odds_ratio_causal/odds_ratio_allestr
