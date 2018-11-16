#!/bin/bash
./snp_str_ld_calculator_MG.py --str-vcf /storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz --snp-vcf /storage/resources/datasets/gtex/59533/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz --loci-file /storage/szfeupe/Runs/650GTEx_estr/gwas/tmp/str_gwas_overlaps_ALL.bed --use-info-start --mincount 3 --usefilter --use-gb > /storage/szfeupe/Runs/650GTEx_estr/gwas/tmp/str_gwas_ld_hg19_ALL.tab

