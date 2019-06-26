#!/bin/bash

# Get gene annotation file
./get_gene_annot.sh > gencode.v7.tab

# Build per-tissue master tables
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

for tissue in $TISSUES
do
    echo $tissue
    ./BuildMasterTable.py \
	--hipref /storage/resources/dbase/human/hg19/hg19.hipstr_reference_withmotif_stranded.bed \
	--tissue $tissue \
	--linreg /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab \
	--mashr-beta /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/posterior_betas.tsv \
	--mashr-se /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/posterior_beta_ses.tsv \
	--anova "/storage/mgymrek/gtex-estrs/revision/anova/${tissue}*.tab" \
	--caviar "/storage/mgymrek/gtex-estrs/revision/caviar/batches/output/${tissue}*.tab" \
	--geneannot gencode.v7.tab \
	--out /storage/mgymrek/gtex-estrs/revision/mastertables/${tissue}_master.tab
done
