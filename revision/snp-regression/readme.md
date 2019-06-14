##### Upload necessary data to AWS: #####

1. Expression datasets for each tissue
```
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
	aws s3 cp /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${tissue}/Corr_Expr.csv s3://gtex-estr/Corr_Expr_${tissue}.csv
done
```

2. Expression annotation file
```
aws s3 cp /storage/resources/dbase/human/hg19/gencode_gene_annotations_hg19.csv s3://gtex-estr/gencode_gene_annotations_hg19.csv
```

3. SNP genotype matrices
```
for chrom in $(seq 1 22)
do
	aws s3 cp /storage/szfeupe/Runs/650GTEx_estr/SNP_Analysis/chr${chrom}.tab s3://gtex-estr/snp_gts_chr${chrom}.tab
done
```

##### Make Dockerfile to run regression analysis (1 job/chrom) #####

```
docker build -t gymreklab/gtex-estrs-snpreg .
```

# Test
```
docker run -it \
       --env AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
       --env AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
       gymreklab/gtex-estrs-snpreg \
       WholeBlood 21
```

##### Set up AWS Batch environment #####
