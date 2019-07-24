#!/bin/bash

# Get overlap of FM-eSTRs with GWAS catalog (for Supp Dataset 3)

./GetGWASOverlap.py \
    --causal /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/ucsc_gwas_catalog_072419_v2.tab \
    --window 50000
