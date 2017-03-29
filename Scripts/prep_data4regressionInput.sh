#!/bin/bash

snp_vcf=/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1/GTEx_Analysis_20150112_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz

#Sample attribute file location

sample_attributes=/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/PhenotypeFiles/phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt.gz

#Define output here
output=GTEx_WGS_SNP_CNV.genotype-calls

#########Extract the principal components. Top 20 by default; this can be cahnge to any x with [ --pca x ] 
	##QC variants 
#echo 'Filtering variants that passed filters'
#vcftools --gzvcf $snp_vcf --remove-filtered-all --remove-indels --recode --stdout|  gzip -c > /storage/szfeupe/data/output_PASS_only.vcf.gz


snp_vcf=/storage/szfeupe/data/output_PASS_only.vcf.gz
	##LD pruning
#echo '\n LD-Pruning\n'
#plink --vcf $snp_vcf --memory 100000 --indep-pairwise 200 100 0.2

	##PCA calculation
#echo '\n PCA calculation'
#plink --vcf $snp_vcf --exclude plink.prune.out --maf 0.05 --hwe 0.0001 --pca header --memory 100000 --out ../data/$output
#echo 'END plink'

#########Get sample attributes and first 3 pcs
	##Sample race, gender and ethnicity
#zgrep "^[^#]" $sample_attributes | awk '{print $2,$4,$6,$7 }' > ../data/SampleAttributes.txt
echo '\n----->Race, Ethnicity and gender extracted from main files...'

	##get first 3 Eigen vectors from  file
ext='.eigenvec'
cut -c -9 ../data/$output$ext > tmp1
awk '{print $3,$4,$5}' ../data/$output$ext >tmp2
echo '' > ../data/Eigenvec.1.2.3.txt
paste tmp1 tmp2 >> ../data/Eigenvec.1.2.3.txt
echo '\n----->Samples first 3 Eigenvec retrieved...'

	## Add race to Eigenvec file
#cleanup
echo '\n----->Cleanup'
rm tmp1
rm tmp2
#rm plink.prune* 
#matching IDs in Eignvec and Attributes: We only take the race($3 of SampleAttribute)

echo 'ID      PC1     PC2     PC3     PCC' > Eigenvec_race.txt
awk 'FNR==NR {a[$1]=$0; next}; $1 in a {print a[$1]"   "$3}' ../data/Eigenvec.1.2.3.txt ../data/SampleAttributes.txt>> Eigenvec_race.txt

##########Plot the first 2 princ comp.
#Call R script for testing the first 2 pc
echo '\n------> Plotting the PCs: Calling plot_pcs.r'
./plot_pcs.r


