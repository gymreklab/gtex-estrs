#!/bin/bash

source params.sh

SAMPLE=$1

HipSTR \
    --bams /storage/gtex-data/wgs/${sample}.bam \
    --fasta /mnt/resources/Homo_sapiens_assembly19.fasta \
    --regions /mnt/resources/GRCh37.hipstr_reference.bed \
    --min-reads 10 \
    --stutter-in /mnt/resources/stutter_logs_0928.txt \
    --str-vcf /storage/vcfs/${sample}.vcf.gz
tabix -p vcf /storage/vcfs/${sample}.vcf.gz
