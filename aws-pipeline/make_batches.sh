#!/bin/bash

source params.sh

rm -f batches/*

split -l ${BATCHSIZE} -d -a 3 gtex_ids.txt ${GTEXBATCHES}/gtex"."

# Sync with S3
aws s3 sync ${GTEXBATCHES} s3://gtex-hipstr/batches/
