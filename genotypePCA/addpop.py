#!/usr/bin/env python

"""
Add population labels

Usage: ./addpop.py indfile popfile > newindfile
"""

import gzip
import sys

try:
    indfile = sys.argv[1]
    sampfile = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Load sample -> population label dictionary
sample_to_pop = {}
labels = ["Asian", "AfricanAmerican", "European", "Amerindian", "NotReported", "Unknown"]
nums = [1, 2, 3, 4, 98, 99]
num_to_pop = dict(zip(nums, labels))
with gzip.open(sampfile, "r") as f:
    for line in f:
        if line.startswith("#"): continue
        if line.strip() == "": continue
        if line.startswith("dbGaP_Subject_ID"): continue
        items = line.strip().split()
        try:
            sample = items[1]
            pop = int(items[5])
            sample_to_pop[sample] = num_to_pop[int(pop)]
        except IndexError: continue

with open(indfile, "r") as f:
    for line in f:
        items = line.strip().split()
        sampid = "-".join(items[0].split("-")[0:2])
        pop = sample_to_pop.get(sampid, "NA")
        items[2] = pop
        sys.stdout.write("\t".join(items)+"\n")
        
