#!/bin/bash

#GTEXDIR=/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/
GTEXDIR=/storage/resources/datasets/gtex/59533/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/
#VCFFILE=${GTEXDIR}/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1/GTEx_Analysis_20150112_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz
VCFFILE=${GTEXDIR}/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
SAMPLEFILE=${GTEXDIR}/PhenotypeFiles/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz
#SAMPLEFILE=/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/PhenotypeFiles/phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt.gz
HMAPFILE=/storage/resources/datasets/hapmap3/hapmap3_b37
MINMAF=0.05

OUTDIR=/storage/mgymrek/gtex/genotypePCA/
PREFIX=${OUTDIR}/GTEx_wgs_650h
PREFIX2=${OUTDIR}/Hapmap3
ALLPREFIX=${OUTDIR}/GTEx_hapmap_merged_650h

# Convert VCF to ped format
vcftools --gzvcf ${VCFFILE} --plink --out ${PREFIX} --remove-indels --remove-filtered-all --maf ${MINMAF}

# Make map files same format
cat ${PREFIX}.map | awk '{print $1 "\t" $1"_"$4 "\t" $3 "\t" $4}' > ${PREFIX}.map.tmp
mv ${PREFIX}.map.tmp ${PREFIX}.map
cp ${HMAPFILE}.ped ${PREFIX2}.ped
cat ${HMAPFILE}.map | awk '{print $1 "\t" $1"_"$4 "\t" $3 "\t" $4}' > ${PREFIX2}.map.tmp
mv ${PREFIX2}.map.tmp ${PREFIX2}.map

# Get only biallelic
plink --file ${PREFIX} --biallelic-only --out ${PREFIX}.biallelic --make-bed --exclude ${ALLPREFIX}-merge.missnp
plink --file ${PREFIX2} --biallelic-only --out ${PREFIX2}.biallelic --make-bed --exclude ${ALLPREFIX}-merge.missnp 

# Merge ped files
plink --bfile ${PREFIX}.biallelic --bmerge ${PREFIX2}.biallelic --make-bed --out ${ALLPREFIX}
# (then go back and remove problematic snps before remerging)

# LD prune
plink --bfile ${ALLPREFIX} --memory 100000 --indep 50 5 2 --out ${ALLPREFIX}
plink --bfile ${ALLPREFIX} --exclude ${ALLPREFIX}.prune.out --maf ${MINMAF} --out ${ALLPREFIX}.pruned \
    --recode \
    --geno 0.05

# Convert to eigstrat format
parfile=convertf_parfile_v2h.txt
echo "genotypename: " ${ALLPREFIX}.pruned.ped > ${parfile}
echo "snpname: " ${ALLPREFIX}.pruned.map >> ${parfile}
echo "indivname: " ${ALLPREFIX}.pruned.ped >> ${parfile}
echo "outputformat: EIGENSTRAT" >> ${parfile}
echo "genotypeoutname: " ${ALLPREFIX}.eigenstratgeno >> ${parfile}
echo "snpoutname: " ${ALLPREFIX}.snp >> ${parfile}
echo "indivoutname: " ${ALLPREFIX}.ind >> ${parfile}
echo "familynames: NO" >> ${parfile}
convertf -p ${parfile}

# Add population group to the ind file
/home/mgymrek/workspace/gtex-estrs/genotypePCA/addpop.py ${ALLPREFIX}.ind ${SAMPLEFILE} > ${ALLPREFIX}.ind.poplabels

# Run smartpca.pl - 1kg only
smartpca.perl \
    -i ${ALLPREFIX}.eigenstratgeno \
    -a ${ALLPREFIX}.snp \
    -b ${ALLPREFIX}.ind.poplabels \
    -k 10 \
    -o ${ALLPREFIX}.1kg.pca \
    -e ${ALLPREFIX}.1kg.evals \
    -p ${ALLPREFIX}.1kg.plot \
    -l ${ALLPREFIX}.1kg.log \
    -w /home/mgymrek/workspace/gtex-estrs/genotypePCA/1kg_pops.txt

# Remove intermediate files
rm ${PREFIX}*
rm ${PREFIX2}*
