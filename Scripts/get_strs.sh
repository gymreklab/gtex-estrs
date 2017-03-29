#!/bin/bash

##USAGE: [./runparallel.sh] or ./get_strs.sh chrNN  where NN= 1 to 22 or Y,X
## --bams [list of bam files]
## --fasta [As I understand this would be the reference seq correcponding to the search domain]
## --region     [This is a requires a BED file containing the STR regions to genotype]
##              BUT What about general WGS screening of STR screening?]
## --stutter-out   [As I understand this will be created along with the output vcf as the model is learned]
## --chrom      [This option is the only available speed option to analyze each chromosome in parallel]
## The file content of "population.list" is the path to each and every WGS included in the calls
##
##Be sure to add the population.list file in the running DIR. The file has in each line the path to an index
ed BAM
##This code is meant to run HipSTR on aligned bam files following the documentation
## STR_DIR is the directory that will contain the output
##Create output directory
#mkdir STR_DIR

############Get population.list
##ls /storage/resources/datasets/gtex/bams/*.bai| sed 's/.bai//' >population.list

############ Fill up the options
chroms=$1
Ref=/storage/szfeupe/refseq/Homo_sapiens_assembly19.fasta       #Reference
out=/storage/szfeupe/Runs/GTEx_estr                             #Output DIR
#REGION=Region_chr1.bed                                         #STR Regions
REGION=/storage/szfeupe/data/Regions_HipSTR_reference.hg19.bed  #STR Regions
str=$out/STR_out.$chroms.vcf.gz                                 #Output filename
filtered=$out/filtered_$chroms.vcf.gz                           #Filtered output

##Run HipSTR###############################

##format bam file to be coma separated
#population=$(sed ':a;N;$!ba;s/\n/,/g' population)

#/home/szfeupe/bin/HipSTR/HipSTR --bam-files population.list --fasta $Ref  --chrom chr1 --regions $REGION --
stutter-out $out/stutter_chr1.txt --str-vcf $str


/usr/bin/HipSTR --bam-files population.list --fasta $Ref  --chrom $chroms --regions $REGION --stutter-out $o
ut/stutter_$chroms --str-vcf $str









