#!/bin/bash

SUPERBATCHPATH=$1

HOMEDIR=/root/

OUTBUCKET=s3://gtex-hipstr/vcfs/

HOMEDIR=/root/

superbatch=$(basename $SUPERBATCHPATH)

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
#    sudo aws ec2 terminate-instances --output table --instance-ids ${INSTANCE_ID}
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

# Download all bam files to EBS storage
sudo mkdir -p /storage/gtex-data
cd /storage/gtex-data/ # go to dbgap directory
for sample in $(cat /storage/tmp/superbatch.txt)
do
    # Download files using aspera
    prefetch --max-size 200G ${sample}
    sam-dump -u sra/${sample}.sra | samtools view -bS > wgs/${sample}.bam
    samtools index wgs/${sample}.bam
    # Run HipSTR
    cd ${HOMEDIR}/source/gtex-estrs/aws-pipeline
    ./run_hipstr_gtex_aws.sh ${sample}
    aws s3 cp /storage/vcfs/${sample}.vcf.gz ${OUTBUCKET}/${sample}_hipstr.vcf.gz
    aws s3 cp /storage/vcfs/${sample}.vcf.gz.tbi ${OUTBUCKET}/${sample}_hipstr.vcf.gz.tbi
    aws s3 cp /storage/vcfs/${sample}.log.txt ${OUTBUCKET}/${sample}_hipstr.log.txt
    # Remove files
    rm -rf sra/${sample}*
    rm -rf wgs/${sample}*
    rm -rf /storage/vcfs/${sample}*
done

terminate

