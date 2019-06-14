#!/bin/bash

CHROM=$1
TISSUE=$2

# Get all the files we need from S3
aws s3 cp s3://gtex-estr/snp_gts_chr${CHROM}.tab /scratch/snp_gts_chr${CHROM}.tab
aws s3 cp s3://gtex-estr/gencode_gene_annotations_hg19.csv /scratch/gencode_gene_annotations_hg19.csv
aws s3 cp s3://gtex-estr/Corr_Expr_${TISSUE}.csv /scratch/Corr_Expr_${TISSUE}.csv

# Run regression analysis
./LinRegAssociationTest.py \
    --chrom ${CHROM} \
    --expr /scratch/Corr_Expr_${TISSUE}.csv \
    --exprannot /scratch/gencode_gene_annotations_hg19.csv \
    --strgt /scratch/snp_gts_chr${CHROM}.tab \
    --distfromgene 100000 \
    --norm \
    --out Lin_Reg_Out_SNPs_${TISSUE}_${CHROM}.tab \
    --tmpdir /tmp
