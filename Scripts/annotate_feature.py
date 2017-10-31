#!/usr/bin/env python
"""
Annotate merged causality score table with distance to TSS and TES of target gene
"""

import numpy as np
import pandas as pd
import sys

try:
    mergefile = sys.argv[1]
    geneannot = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

data = pd.read_csv(mergefile, sep="\t")
annot = pd.read_csv(geneannot)
annot["gene"] = annot["gene.id"]
tmpdata = pd.merge(data, annot[["gene","gene.strand","gene.start","gene.stop"]], on=["gene"])

# Annotate TSS/TES
def GetDistTss(x):
    if x["gene.strand"] == "+":
        return x["best.str.start"] - x["gene.start"]
    elif x["gene.strand"] == "-":
        return -1*(x["best.str.start"] - x["gene.stop"])
    else:
        return float("nan")

def GetDistTes(x):
    if x["gene.strand"] == "+":
        return x["best.str.start"] - x["gene.stop"]
    elif x["gene.strand"] == "-":
        return -1*(x["best.str.start"] - x["gene.start"])
    else:
        return float("nan")

tmpdata["dist.to.tss"] = tmpdata.apply(lambda x: GetDistTss(x), 1)
tmpdata["dist.to.tes"] = tmpdata.apply(lambda x: GetDistTes(x), 1)
tmpdata[["gene","chrom","best.str.start","best.score","best.tissue","best.q","dist.to.tss","dist.to.tes"]].to_csv("GTEx_merged_causality_tsstes_", sep="\t", index=False)

def GetBootstrapCI(data, func, perc):
    numiter = 1000
    vals = []
    for i in range(numiter):
        x = np.random.choice(data, size=data.shape[0], replace=True)
        vals.append(func(x))
    percs = map(lambda x: int(x*100), [(1-perc)/2, perc + (1-perc)/2])
    return np.percentile(vals, q=percs)

# Print summary to stderr
sys.stderr.write("\t".join(["bin.start","bin.end","bin.size", \
                            "perc.tss", "perc.tss", \
                            "enrich.tss", "enrich.tss.low", "enrich.tss.high", \
                            "enrich.tes", "enrich.tes.low", "enrich.tes.high"])+"\n")
MINDIST = 1000
bins = np.percentile(tmpdata["best.score"], q = [0, 50] + list(np.arange(60, 101, 10)))
for i in range(len(bins)-1):
    lb = bins[i]
    ub = bins[i+1]
    x = tmpdata[(tmpdata["best.score"]>lb) & (tmpdata["best.score"]<= ub)]
    feature_tss = x["dist.to.tss"].apply(lambda x: abs(x)<MINDIST)
    feature_tes = x["dist.to.tes"].apply(lambda x: abs(x)<MINDIST)
    count_tss = np.mean(feature_tss)
    count_tes = np.mean(feature_tes)
    low_tss, high_tss = GetBootstrapCI(feature_tss, np.mean, 0.95)
    low_tes, high_tes = GetBootstrapCI(feature_tes, np.mean, 0.95)
    if i == 0:
        base_tss = count_tss
        base_tes = count_tes
    sys.stderr.write("\t".join(map(str, [lb, ub, x.shape[0], \
                                         count_tss, count_tes, \
                                         count_tss*1.0/base_tss, low_tss*1.0/base_tss, high_tss*1.0/base_tss, \
                                         count_tes*1.0/base_tes, low_tes*1.0/base_tes, high_tes*1.0/base_tes]))+"\n")

