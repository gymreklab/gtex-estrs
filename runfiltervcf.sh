#!/bin/bash

OUT=Filtered_$1
IN=/storage/resources/datasets/gtex/vcfs_650/

python Scripts/filter_vcf.py  --vcf $IN$1 --min-call-qual 0.9 --max-call-flank-indel  0.15 --max-call-stutter  0.15 >/storage/szfeupe/Runs/650GTEx_estr/VCFs/$OUT

