#!/bin/bash

# Run mfold on a couple loci - CSTB example
./run_mfold_locus.sh 21:45196326-45196326 

# Run on some causal examples
# Replace this with something cleaner 
# e.g. getting directly from master tables
cat ~/workspace/gtex-estrs/Scripts/Figures/allcestrs3_tmp.tab | \
    grep -v chrom | awk '{print $1 "\t" $2 "\t" $2+1}' | \
    intersectBed -a stdin -b /storage/resources/dbase/human/hg19/hg19.hipstr_reference.bed -wa -wb | \
    awk '($2==$5)' | awk '{print $1":"$2"-"$6}' | sort | uniq | sed 's/chr//' > causal_example_loci.bed

cat causal_example_loci.bed | xargs -n 1 -P 1 ./run_mfold_locus.sh


