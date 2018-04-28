#!/usr/bin/env python

"""
Download BED files for TFBS from ENCODE

Usage:
./download_tfbps.py <accs> <outdir>

e.g. ./download_tfbps.py encode_tfhistone_files.txt /storage/mgymrek/gtex/tfbs/encodedata
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

h = 1
with open(accfile, "r") as f:
    for line in f:
        if h:
            h = 0
            continue
        fileurl = line.strip()
        expurl = fileurl.split("@")[0]
        response = requests.get(expurl, headers=HEADERS)
        response_json_dict = response.json()
        acc = response_json_dict["dataset"].split("/")[2]
        URL = "https://www.encodeproject.org/experiments/%s/"%acc
        response = requests.get(URL, headers=HEADERS)
        response_json_dict = response.json()
        if "target" not in response_json_dict.keys():
            sys.stderr.write("Couldn't find info for %s\n"%acc)
            continue
        target = response_json_dict["target"]["title"].split()[0]
        celltype = response_json_dict["biosample_summary"]
        if celltype == "K562": continue
        output_file = os.path.join(outdir, "%s_%s.bed.gz"%(target, celltype))
        cmd = "wget -O %s %s"%(output_file, fileurl)
        os.system(cmd)
