#!/bin/bash

# Index TFBS regions with Giggle
./index_tfs.sh

# Get overlaps for each tissue
./tf_overlaps.sh WholeBlood
