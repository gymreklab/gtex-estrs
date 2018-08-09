#!/bin/bash

source params.sh

GWASBED=$1 # chrom, start, rsid, locusname, pval
PREFIX=$2

##### Get all STRs within 50kb of a GWAS catalog SNP
echo "Get all nearby STRs...."
cat ${GWASBED} | grep -v chrom | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' > ${OUTDIR}/tmp/${PREFIX}.bed
cat ${HIPREF} | sed 's/chr//' | awk -v"window=$WINDOW" '{print $1 "\t" ($2<window?0:$2-window) "\t" $3+window "\t" $0}' | \
    intersectBed -a stdin -b ${OUTDIR}/tmp/${PREFIX}.bed -wa -wb | cut -f 1-3 --complement  > ${OUTDIR}/tmp/str_gwas_overlap_${PREFIX}.bed
cat ${OUTDIR}/tmp/str_gwas_overlap_${PREFIX}.bed | cut -f 1,2,7 | uniq > ${OUTDIR}/tmp/str_gwas_coords_${PREFIX}.tab

##### Get SNP-STR LD
echo "Get SNP-STR LD..."
~/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
    --str-vcf ${STRVCF} \
    --snp-vcf ${SNPVCF} \
    --loci-file ${OUTDIR}/tmp/str_gwas_coords_${PREFIX}.tab \
    --use-info-start --mincount 3 --usefilter --use-gb > ${OUTDIR}/tmp/str_gwas_ld_hg19_${PREFIX}.tab

##### Combine
echo "Combine..."
cat ${OUTDIR}/tmp/str_gwas_ld_hg19_${PREFIX}.tab | grep -v locus1 | sed 's/:/\t/g' | \
    awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' | intersectBed -a stdin -b ${OUTDIR}/tmp/str_gwas_overlap_${PREFIX}.bed -wa -wb | \
    awk '(($2==$15) && ($7==$20))' | cut -f 4,5,12,17,18,20,24- > ${OUTDIR}/str_gwas_ld_COMBINED_${PREFIX}.tab

##### Combine with eSTRs
cat ${OUTDIR}/str_gwas_ld_COMBINED_${PREFIX}.tab | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' | sed 's/chr//' > ${OUTDIR}/tmp/del1${PREFIX}
cat ${CAUSAL} | grep -v gene | awk '{print $2 "\t" $3 "\t" $3+1 "\t" $0}' | sed 's/chr//' > ${OUTDIR}/tmp/del2${PREFIX}
echo "chrom,start,LD,period,motif,snp.start,rsid,locus,pval,gene,best.score,bes.q,best.variant,best.tissue,e.tissues" | sed 's/,/\t/g' \
    > ${OUTDIR}/str_gwas_ld_COMBINED_eSTR_${PREFIX}.tab 
intersectBed -a ${OUTDIR}/tmp/del1${PREFIX} -b ${OUTDIR}/tmp/del2${PREFIX} -wa -wb  | cut -f 4-12,16,19- | uniq \
    >> ${OUTDIR}/str_gwas_ld_COMBINED_eSTR_${PREFIX}.tab 

##### Get candidates
echo "Get candidates..."
cat ${OUTDIR}/str_gwas_ld_COMBINED_eSTR_${PREFIX}.tab | \
    awk -F"\t" '($3>0.2 && $12<0.1 && $11>0)' | grep -v "SNP_" | sort -k 3 -g -r
