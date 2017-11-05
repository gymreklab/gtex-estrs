#!/usr/bin/env python

"""
Recalculate and merge causality scores from all tissues
Usage: ./merge_causality_scores.py <basedir> <tissuelist> <scoretype>
Scoretype can be:
- causality
- posterior
"""

import os
import pandas as pd
import sys

SCORETYPES = ["causality", "posterior"]
try:
    basedir = sys.argv[1]
    tissues = sys.argv[2].split(",")
    scoretype = sys.argv[3]
    if scoretype not in SCORETYPES:
        sys.stderr.write("__doc__")
        sys.exit(1)
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Load data for each tissue, keep list of genes
tissue_data = {}
genes = set()
for t in tissues:
    data = pd.read_csv(os.path.join(basedir, t, "Master.table"), sep="\t")
    # Recalculate causality score
    if scoretype == "causality":
        data["cis_str_h2"] = data["cis_str_h2"].apply(lambda x: min(1, max(x, 10**-6)))
        data["cis_snp_h2"] = data["cis_snp_h2"].apply(lambda x: min(1, max(x, 10**-6)))
        data["cis_h2"] = data["cis_str_h2"] + data["cis_snp_h2"]
        data["score"] = data.apply(lambda x: x["best.str.score"]*(x["cis_str_h2"]/(x["cis_h2"])), 1)
    elif scoretype == "posterior":
        data["score"] = data["best.str.score"]
    else: data["score"] = -1
    tissue_data[t] = data[["gene","chrom","best.str.start","score","qvalue"]]
    genes = genes.union(set(data["gene"]))

genes = list(genes)
d_chrom = []
d_pos = []
d_score = []
d_tissue = []
d_qval = []

for gene in genes:
    best_score = -1
    best_tissue = "NA"
    best_str = "NA"
    best_q = "NA"
    chrom = "NA"
    for t in tissues:
        x = tissue_data[t]
        x = x[x["gene"]==gene]
        if x.shape[0] == 1:
            score = x["score"].values[0]
            chrom = x["chrom"].values[0]
            start = x["best.str.start"].values[0]
            q = x["qvalue"].values[0]
            if score > best_score:
                best_score = score
                best_tissue = t
                best_str = start
                best_q = q
                chrom = chrom
    d_chrom.append(chrom)
    d_pos.append(best_str)
    d_score.append(best_score)
    d_tissue.append(best_tissue)
    d_qval.append(best_q)

df = pd.DataFrame({"gene": genes,
                   "chrom": d_chrom,
                   "best.str.start": d_pos,
                   "best.score": d_score,
                   "best.tissue": d_tissue,
                   "best.q": d_qval})
df[["gene","chrom","best.str.start","best.score","best.tissue","best.q"]].to_csv(sys.stdout, sep="\t", index=False)