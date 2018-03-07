#!/bin/bash

source params.sh

# Sort input files
#mkdir -p ${TFDIRSORT}
#/home/mgymrek/workspace/giggle/scripts/sort_bed "${TFDIR}/*.bed.gz" ${TFDIRSORT} 4

giggle index -i "${TFDIRSORT}/*.gz" -s -f -o ${INDEXDIR}
