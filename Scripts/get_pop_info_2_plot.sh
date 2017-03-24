#!/bin/bash

##retrieve population info

#Sample attribute file location
sample_attributes=/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/PhenotypeFiles/phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt.gz

grep '^#' $sample_attrib | awk '{print $2,$4,$6,$7 }' > ../data/SampleAttributes.txt
echo 'Race, Ethnicity and gender extracted from main files...'
##get first 3 Eigen vectors from  file 

cut -c -9 ../SampleAttributes.txt > ../data/Eigen1.2.3.txt
echo 'Sample attribute IDs retrieved...'

paste Eigen1.2.3.txt < (awk '{print $3,$4,$5}' ../data/GTEx_WGS_SNP_CNV.genotype-calls_chr1.eigenvec )
echo 'IDs and first 3 Eigen vectors combined.... Done'
##
##
