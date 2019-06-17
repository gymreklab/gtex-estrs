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
docker push gymreklab/gtex-estrs-snpreg
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

```
aws batch create-compute-environment \
    --compute-environment-name c42xlarge \
    --type MANAGED \
    --state ENABLED \
    --compute-resources file://batch-small.json \
    --service-role arn:aws:iam::369425333806:role/service-role/AWSBatchServiceRole

aws batch create-job-queue \
    --job-queue-name gtex-c42xlarge \
    --state ENABLED \
    --priority 100 \
    --compute-environment-order order=1,computeEnvironment=c42xlarge

aws batch register-job-definition \
    --job-definition-name gtex-snpreg-job \
    --type container \
    --container-properties file://gtex-snpreg-job.json
```

##### Run AWS jobs #####

```
tissue="Nerve-Tibial"
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 #$(seq 1 20) # 21 22
do
cmd="aws batch submit-job \
    --job-name NerveTibial-${chrom} \
    --job-queue gtex-c42xlarge \
    --job-definition gtex-snpreg-job:6 \
    --container-overrides 'command=[\"${tissue}\",\"${chrom}\"]'"
sh -c "${cmd}"
done 
```