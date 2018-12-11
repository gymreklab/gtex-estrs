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

# Combine to one output file
cat *.energies.tab | sort -k 1,1 -k2,2n > causal_energies.tab

# Overlap with score, beta, and motif
cat ~/workspace/gtex-estrs/Scripts/Figures/allcestrs3_tmp.tab | grep -v chrom | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $4 "\t" $6 "\t" $7}' | intersectBed -a stdin -b /storage/resources/dbase/human/hg19/hg19.hipstr_reference.bed -wa -wb | awk '($2==$8)' | awk '{print $1"_"$2"-"$9 "\t" $4 "\t" $5 "\t" $6}'  | sort -k 1,1 | sed 's/chr//' > causal_loci_info.tab
cat causal_energies.tab | datamash -g 1 ppearson 2:4 min 4 max 4| sort -k 1,1 > causal_energies_corr.tab

# Now for control set
cat ~/workspace/gtex-estrs/Scripts/Figures/allstrs_tmp.csv | sed 's/,/\t/g' | awk '{print $3 "\t" $4 "\t" $4+1}' | grep -v chrom | \
    intersectBed -a stdin -b /storage/resources/dbase/human/hg19/hg19.hipstr_reference.bed -wa -wb | \
    awk '($2==$5)' | awk '{print $1":"$2"-"$6}' | sort | uniq | sed 's/chr//' > all_example_loci.bed
cat all_example_loci.bed | xargs -n 1 -P 1 ./run_mfold_locus.sh
