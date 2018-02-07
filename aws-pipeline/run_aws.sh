#!/bin/bash

source params.sh

AWS_ACCESS_KEY=$(cat ~/.aws/credentials | grep "aws_access" | cut -f 3 -d' ' | head -n 1)
AWS_SECRET_KEY=$(cat ~/.aws/credentials | grep "aws_secret" | cut -f 3 -d' ' | head -n 1)

for batch in $(ls -l ${GTEXBATCHES} | grep -v total | awk '{print $NF}')
do
    batchpath=s3://gtex-hipstr/batches/${batch}
    ./launch_aws.sh ${batchpath} ${AWS_ACCESS_KEY} ${AWS_SECRET_KEY} micro_key
done
