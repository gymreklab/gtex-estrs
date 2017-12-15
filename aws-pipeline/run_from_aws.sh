#!/bin/bash

SUPERBATCHPATH=$1

HOMEDIR=/root/

OUTBUCKET=s3://gtex-hipstr/vcfs/

superbatch=$(basename $SUPERBATCHPATH)

HOMEDIR=/root/

usage()
{
    BASE=$(basename -- "$0")
    echo "Run HipSTR on GTEx
Usage:
    $BASE <superbatch>

Does the following:
1. Set up
2. Download necessary files
3. Create jobs
4. Run those jobs
5. Upload results to S3 bucket
6. Terminate
"
    terminate
    exit 1
}

terminate() {
    INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)
    # Get log
    sudo aws s3 cp --output table /var/log/cloud-init-output.log ${OUTBUCKET}/log/${superbatch}.log
    # Terminate instance
    echo "Terminating instance ${INSTANCE_ID}"
    sudo aws ec2 terminate-instances --output table --instance-ids ${INSTANCE_ID}
    exit 1 # shouldn't happen
}

test -z ${SUPERBATCHPATH} && usage

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    terminate
    exit 1
}

# Set of directories for inputs/outsputs on mounted EBS drive
sudo mkfs -t ext4 /dev/xvdf
sudo mkdir /storage
sudo mount /dev/xvdf /storage/
sudo chmod 777 /storage/

# Download files
sudo mkdir -p /storage/tmp || die "Could not make tmp directory"
aws s3 cp ${SUPERBATCHPATH} /storage/tmp/superbatch.txt

# Make directory for inputs/outputs
sudo mkdir -p /storage/vcfs || die "Could not make vcfs directory"

# Get github for gtex project
cd ${HOMEDIR}/source
git clone https://github.com/gymreklab/gtex-estrs

# Run each job
source ${HOMEDIR}/source/gtex-estrs/aws-pipeline/params.sh
for sample in $(cat /storage/tmp/superbatch.txt)
do
    cd /storage/gtex-data/ # go to dbgap directory
    # Download files using aspera
    prefetch --max-size 200G ${sample}
    sam-dump -u sra/${sample}.sra | samtools view -bS > wgs/${sample}.bam
    samtools index wgs/${sample}.bam
    # Run HipSTR
    cd ${HOMEDIR}/source/gtex-estrs/aws-pipeline
    ./run_hipstr_gtex_aws.sh ${sample} 
    # Upload result to S3
    aws s3 cp /storage/vcfs/${sample}_hipstr.vcf.gz ${OUTBUCKET}/${sample}_hipstr.vcf.gz
    aws s3 cp /storage/vcfs/${sample}_hipstr.vcf.gz.tbi ${OUTBUCKET}/${sample}_hipstr.vcf.gz.tbi
    # Delete result and bam
    rm -f wgs/${sample}*
    rm -f /storage/vcfs/${sample}*
done

terminate

