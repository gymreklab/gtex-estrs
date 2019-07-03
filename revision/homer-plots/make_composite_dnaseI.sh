#!/bin/bash

# TODO not done
STRSETS=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/composite
TAGDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/tagdirs
DNASE="${TAGDIR}/DNAse_GM12878 ${TAGDIR}/dnase-heart-ENCFF702IJE ${TAGDIR}/dnase-lung-ENCFF803ZER ${TAGDIR}/dnase-seq-adultheart-ENCFF185ISK ${TAGDIR}/dnase-tibial_artery-ENCFF223EKZ ${TAGDIR}/dnase-tibial_nerve-ENCFF226ZCG"

for period in $(seq 1 6)
do
    /storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/ALLCAUSAL_period${period}.bed \
	hg19 -size 1000 -hist 1 -d ${DNASE} > ${OUTDIR}/ALLCAUSAL_dnase_period${period}.bed
    /storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/ALLSTRs_period${period}.bed \
	hg19 -size 1000 -hist 1 -d ${DNASE} > ${OUTDIR}/ALLSTRs_dnase_period${period}.bed
    for tissue in Heart-LeftVentricle Lung Artery-Tibial Nerve-Tibial
    do
	/storage/resources/source/homer/bin/annotatePeaks.pl \
	    ${STRSETS}/${tissue}_period${period}.bed \
	    hg19 -size 1000 -hist 1 -d ${DNASE} > ${OUTDIR}/${tissue}_dnase_period${period}.bed
    done
done
