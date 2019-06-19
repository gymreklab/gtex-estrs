#!/usr/bin/env python3

"""
Find:
max Z per unit
max Z per gene

Output remaining units with |Z|>=thresh
"""

import os
import pandas as pd
import sys

workdir = sys.argv[1]
prefix = sys.argv[2]

ZTHRESH = 4

# Read in zvals
zvals = pd.read_csv(os.path.join(workdir, "output-%s"%prefix, "zscores.tsv"), sep="\t", index_col=0)

# Compute unit-level and gene-level max
zmax = pd.DataFrame({"ID": zvals.index, "zmax": abs(zvals.apply(max, 1)), "gene": [item.split("_")[0] for item in zvals.index]})
genemax = zmax.groupby("gene", as_index=False).agg({"zmax":max})
genemax.columns = ["gene","genemax"]
zmax = pd.merge(zmax, genemax, on=["gene"])

# Output significant hits, and subset of z-scores for sig hits
sig = zmax[zmax.apply(lambda x: x["zmax"]==x["genemax"] and x["genemax"]>=ZTHRESH, 1)]
sig[["ID","zmax"]].to_csv(os.path.join(workdir, "output-%s"%prefix, "sig.tsv"), sep="\t", index=False)

zsig = zvals.loc[sig["ID"]]
zsig.to_csv(os.path.join(workdir, "output-%s"%prefix, "zscores-sigonly.tsv"), sep="\t", index=True)

# Summarize
print("Sig genes: %s"%sig.shape[0])
print("Sig units: %s"%(zmax[zmax["zmax"]>=ZTHRESH].shape[0]))

