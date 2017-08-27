#!/bin/bash

set -e

# Merge causality scores
BASEDIR=/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/
TISSUES=Adipose-Subcutaneous,Artery-Tibial,Cells-Transformedfibroblasts,Esophagus-Mucosa,Lung,Muscle-Skeletal,WholeBlood
OUTDIR=/storage/mgymrek/gtex/causality
GENEANNOT=/storage/szfeupe/Runs/GTEx_estr/Gene_Exp_Annotation.txt

./merge_causality_scores.py $BASEDIR $TISSUES causality > ${OUTDIR}/GTEx_merged_causality.tab
./merge_causality_scores.py $BASEDIR $TISSUES posterior > ${OUTDIR}/GTEx_merged_causality_post.tab

# Make feature tables (TSS/TES)
./annotate_feature_tsstes.py \
    ${OUTDIR}/GTEx_merged_causality.tab \
    ${GENEANNOT} > ${OUTDIR}/features/GTEx_merged_causality_tsstes.tab

./annotate_feature_tsstes.py \
    ${OUTDIR}/GTEx_merged_causality_post.tab \
    ${GENEANNOT} > ${OUTDIR}/features/GTEx_merged_causality_tsstes_post.tab
