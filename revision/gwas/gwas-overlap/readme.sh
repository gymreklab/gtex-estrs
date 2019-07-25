#!/bin/bash

# Get overlap of FM-eSTRs with GWAS catalog (for Supp Dataset 3)
./GetGWASOverlap.py \
    --estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/ucsc_gwas_catalog_072419_v2.tab \
    --window 50000 \
    --caviar 0.3 \
    --prefix NCBI-072419 \
    --outdir /storage/mgymrek/gtex-estrs/revision/figures/ \
    --str /storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz \
    --snp /storage/resources/datasets/gtex/59533/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz

# See ./prepare_gwas_hits.sh and ./get_gwas_candidates.sh for trait-specific