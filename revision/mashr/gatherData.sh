#!/bin/bash

testrun=false
testrun2Chroms=true

#example where I learned this from
# ls -ld /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/* | grep '^d' | sed -e 's/\s\+/ /g' | cut -f 9 -d ' '

#TODO: check this!
tissuetypes='Adipose-Subcutaneous
Adipose-Visceral(Omentum)
Artery-Aorta
Artery-Tibial
Brain-Caudate(basalganglia)
Brain-Cerebellum
Cells-Transformedfibroblasts
Esophagus-Mucosa
Esophagus-Muscularis
Heart-LeftVentricle
Lung
Muscle-Skeletal
Nerve-Tibial
Skin-NotSunExposed(Suprapubic)
Skin-SunExposed(Lowerleg)
Thyroid
WholeBlood
'

#Don't seem to have tables
#'Kidney-Cortex Liver '


#escape the tissue types
tissuetypes=$(printf '%q ' $tissuetypes)

if [ "$testrun" = true ]; then
	#only use three tissues - two similar, one not, in testing
	testtissues='Artery-Aorta Artery-Tibial Lung'
	tissuetypes=$testtissues
	workdir='testrun'
elif [ "$testrun2Chroms" = true ]; then
	#only use three tissues - two similar, one not, in testing
	testtissues='Artery-Aorta Artery-Tibial Lung'
	tissuetypes=$testtissues
	workdir='testrun2Chroms'
else
	workdir='fullrun'
fi

mkdir -p ${workdir}/output
mkdir -p ${workdir}/input
mkdir -p ${workdir}/intermediate

#Example command
#cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/WholeBlood/Master.table | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $27}' | head
#relevant columns are 7,9,10,26,27
tissuedir=/storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue
for tissue in $tissuetypes ; do
	command="cat ${tissuedir}/${tissue}/Master.table"' | awk -F"\t" '"'"'{print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $27}'"'"

	if [ "$testrun" = true ]; then
		#only use chromosome 21 in testing
		command="${command} | grep -P '(chr21\t)|(chrom)'"
	elif [ "$testrun2Chroms" = true ]; then
		command="${command} | grep -P '(chr21\t)|(chrom)|(chr22\t)'"
	fi
	command="${command} > ./${workdir}/input/$tissue.table"
	eval $command
done


