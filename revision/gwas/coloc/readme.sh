# See make_coloc_summstats.sh to get summ stats in necessary format
# SCZPGC HeightYengo IBDHuang IntelligenceSavageJensen AlzheimersIGAP
OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/results/

./run_trait.sh SCZPGC cc 0.33
./run_trait.sh HeightYengo quant 695647
./run_trait.sh IBDHuang cc 0.33
./run_trait.sh IntelligenceSavageJensen quant 269720
./run_trait.sh AlzheimersIGAP cc 0.33

for trait in SCZPGC HeightYengo
do
    ./summarize_coloc.sh ${trait} > ${OUTDIR}/${trait}_coloc_results.tab
done

#    cat ${OUTDIR}/${trait}_coloc_results.tab | grep -v gene | awk '($12>0.5)' | \
#	sort -k6,6g | awk '{print $1":"$4 "\t" $2 "\t" $3 "\t" $5 "\t" $6"\t" $7 "\t" $12}'
