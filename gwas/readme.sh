#!/bin/bash

# Get GWAS overlaps
./get_gwas_overlaps.sh

# Filter gwas tables and get supp dataset 2
./filter_gwas.sh GWASCAT2

