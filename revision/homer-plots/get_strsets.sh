#!/bin/bash

OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
MASTER=/storage/mgymrek/gtex-estrs/revision/mastertables/
CAUSAL=/storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv
SCORE=0.3 # caviar threshold

# Get all STRs by period
for period in $(seq 1 6)
do
    cat ${MASTER}/* | \
	grep -v chrom | awk -v"period=$period" '(length($9)==period) {print $1 "\t" $3 "\t" $8}'| \
	sort | uniq > ${OUTDIR}/ALLSTRs_period${period}.bed
done
exit 1
# Get aggregate eSTRs by period
for period in $(seq 1 6)
do
    cat ${CAUSAL} | csvcut -t -c chrom,str.start,str.end,str.motif.forward,score | \
	csvformat -T | grep -v chrom | \
	awk -v"period=$period" -v"score=$SCORE" \
	'(length($4)==period && $5>=score) {print $1 "\t" $2 "\t" $3}' | \
	sort | uniq > ${OUTDIR}/ALLCAUSAL_period${period}.bed
done

# Get per-tissue eSTRs by period
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
    for period in $(seq 1 6)
    do
	cat ${MASTER}/${tissue}_master.tab | \
	    csvcut -t -c chrom,str.start,str.end,str.motif.forward,caviar.str.score | \
	    csvformat -T | grep -v chrom | \
	    awk -v"period=$period" -v"score=$SCORE" \
	    '(length($4)==period && $5>=score) {print $1 "\t" $2 "\t" $3}' | \
	    sort | uniq > ${OUTDIR}/${tissue}_period${period}.bed
    done
done
