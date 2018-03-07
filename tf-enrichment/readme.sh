#!/bin/bash

source params.sh

# Index TFBS regions with Giggle
#./index_tfs.sh

for t in $TISSUES
do
    echo "Finding overlap for tissue $t"
    # Get overlaps for each tissue - TFBS 
    ./get_overlaps.sh ${t} ${TFINDEX} tfbs
    # Get overlaps for each tissue - EpigenomeRoadmap
    ./get_overlaps.sh ${t} ${ROADMAPINDEX} roadmap
done
