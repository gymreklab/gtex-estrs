#!/bin/bash

source params.sh

# Index TFBS +RNABP regions with Giggle
#./index_tfs.sh
#./index_rnabps.sh

for t in $TISSUES
do
    echo "Finding overlap for tissue $t"
    # Get overlaps for each tissue - TFBS 
#    ./get_overlaps.sh ${t} ${TFINDEX} tfbs
    # Get overlaps for each tissue - EpigenomeRoadmap
#    ./get_overlaps.sh ${t} ${ROADMAPINDEX} roadmap
    # RNABP
    ./get_overlaps.sh ${t} ${RNABPINDEX} rnabp
done
