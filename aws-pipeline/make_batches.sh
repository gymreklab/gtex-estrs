#!/bin/bash

source params.sh

rm -f ${GTEXBATCHES}/*

split -l ${BATCHSIZE} -d -a 3 gtex_ids_todo.txt ${GTEXBATCHES}/gtex"."

# Sync with S3
aws s3 sync ${GTEXBATCHES} s3://gtex-hipstr/batches/
