#!/bin/bash

aws ec2 describe-instances --output text --filters "Name=instance-type,Values=c4.xlarge" "Name=instance-state-name,Values=running" | grep "^INSTANCES" | cut -f 16 > instances.txt

rm -f gtex_aws_progress.txt
touch gtex_aws_progress.txt
while IFS='' read -r line || [[ -n "$line" ]]; do
    instance=$line
    scp -o  StrictHostKeyChecking=no -i ~/keys/micro_key.pem ubuntu@$instance:/var/log/cloud-init-output.log tmp.txt
    batch=$(cat tmp.txt | grep "superbatch" | cut -d' ' -f 9 | awk -F"/" '{print $NF}')
    num=$(cat tmp.txt | grep "checking" | wc -l)
    echo $instance $batch $num >> gtex_aws_progress.txt
done < "instances.txt"

aws s3 ls s3://gtex-hipstr/vcfs/ | grep ".vcf.gz$" | wc -l
