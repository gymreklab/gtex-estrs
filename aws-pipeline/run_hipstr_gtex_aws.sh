#!/bin/bash

SAMPLE=$1

HipSTR \
    --bams /storage/gtex-data/wgs/${SAMPLE}.bam \
    --fasta /mnt/resources/Homo_sapiens_assembly19.fasta \
    --regions /mnt/resources/GRCh37.hipstr_reference.bed \
    --min-reads 5 \
    --str-vcf /storage/vcfs/${SAMPLE}.vcf.gz \
    --stutter-in /mnt/resources/stutter_logs_0928.txt \
    --log /storage/vcfs/${SAMPLE}.log.txt
tabix -p vcf /storage/vcfs/${SAMPLE}.vcf.gz
