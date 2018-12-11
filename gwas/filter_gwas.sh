#!/bin/bash

source params.sh

PREFIX=$1

cat ${OUTDIR}/str_gwas_ld_COMBINED_eSTR_${PREFIX}.tab | \
    awk -v"mincausal=$MINCAUSAL" -F"\t" '($10>mincausal)' | \
    sort -k1,1 -k2,2n -k 7,7 | \
    awk -F "\t" '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' | \
    uniq | grep -v beta | \
    datamash --header-in -g 1,2,3,4,5,6 unique 7 unique 8 max 9 > \
    ${OUTDIR}/${PREFIX}_candidates.tab

# How many unique FM-eSTRs
echo "unique FM-eSTRS " $(cat ${OUTDIR}/${PREFIX}_candidates.tab | cut -f 1,2 | sort | uniq | wc -l)

# How many in moderate, strong LD
echo "moderate LD " $(cat ${OUTDIR}/${PREFIX}_candidates.tab | awk '($3>0.1)' | cut -f 1,2 | sort | uniq | wc -l)
echo "strong LD " $(cat ${OUTDIR}/${PREFIX}_candidates.tab | awk '($3>0.8)' | cut -f 1,2 | sort | uniq | wc -l)

# How many is the lead SNP in the STR itself
cat ${OUTDIR}/${PREFIX}_candidates.tab | awk '{print "chr"$1 "\t" $2 "\t" $2+1 "\t" $0}' | \
    intersectBed -a $HIPREF -b stdin -wa -wb | \
    cut -f 4-10 --complement | \
    awk '($6>=($2-1) && $6<=($3+1))' > ${OUTDIR}/${PREFIX}_candidates_exact.tab
echo "Exact hit " $(cat ${OUTDIR}/${PREFIX}_candidates_exact.tab | cut -f 1,2 | sort | uniq | wc -l)
