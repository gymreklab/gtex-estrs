#!/usr/bin/env python

"""
Download BED files for RNABP sites from clipseq

Usage:
./download_eclip_rnabs.py <accfile> <outdir>

e.g. ./download_eclip_rnabs.py eclip_accessions.txt /storage/mgymrek/gtex/rnabps/encodedata
"""

import os
import requests, json
import sys

try:
    accfile = sys.argv[1]
    outdir = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

with open(accfile, "r") as f:
    for line in f:
        if line.startswith("RBP"): continue
        acc = line.strip().split()[2]
        URL = "https://www.encodeproject.org/experiments/%s/"%acc
        response = requests.get(URL, headers=HEADERS)
        response_json_dict = response.json()
        if "target" not in response_json_dict.keys():
            sys.stderr.write("Couldn't find info for %s\n"%acc)
            continue
        target = response_json_dict["target"]["title"].split()[0]
        celltype = response_json_dict["biosample_summary"]
        files = response_json_dict["files"]
        for fdata in files:
            if fdata["file_type"] == "bed narrowPeak" and \
               fdata["assembly"] == "hg19" and \
               fdata["biological_replicates"][0] == 1:
                download_file = "https://www.encodeproject.org" + fdata["href"]
                output_file = os.path.join(outdir, "%s_%s.bed.gz"%(target, celltype))
                cmd = "wget -O %s %s "%(output_file, download_file)
                print target, celltype
                os.system(cmd)

