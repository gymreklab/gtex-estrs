#!/bin/bash

set -e

SAMPLE=$1

#    --stutter-in /mnt/resources/stutter_logs_0928.txt \

HipSTR \
    --bams /storage/gtex-data/wgs/${SAMPLE}.bam \
    --fasta /storage/resources/Homo_sapiens_assembly19.fasta \
    --regions /storage/resources/GRCh37.hipstr_reference.bed \
    --min-reads 5 \
    --str-vcf /storage/vcfs/${SAMPLE}.vcf.gz \
    --def-stutter-model \
    --log /storage/vcfs/${SAMPLE}.log.txt
tabix -p vcf /storage/vcfs/${SAMPLE}.vcf.gz

exit 0
