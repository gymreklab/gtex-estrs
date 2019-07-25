#!/bin/bash

trait=$1
type=$2
n=$3

# Get coloc files
./get_coloc_trait.sh ${trait} ${type} ${n}
