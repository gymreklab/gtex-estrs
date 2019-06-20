#!/usr/bin/env python3

import os
import pandas as pd
import sys

workdir = sys.argv[1]
prefix = sys.argv[2]

ZTHRESH = 3

# Read in zvals
zvals = pd.read_csv(os.path.join(workdir, "output-%s"%prefix, "zscores.tsv"), sep="\t", index_col=0)

# Compute unit-level and gene-level max
zmax = pd.DataFrame({"ID": zvals.index, "zmax": abs(zvals.apply(max, 1)), "gene": [item.split("_")[0] for item in zvals.index]})
genemax = zmax.groupby("gene", as_index=False).agg({"zmax":max})
genemax.columns = ["gene","genemax"]
zmax = pd.merge(zmax, genemax, on=["gene"])

# Get significant gene-level hits
sig = zmax[zmax.apply(lambda x: x["zmax"]==x["genemax"] and x["genemax"]>=ZTHRESH, 1)]
sig[["ID","zmax"]].to_csv(os.path.join(workdir, "output-%s"%prefix, "sig.tsv"), sep="\t", index=False)

# Summarize - per tissue
all_estrs = set()
all_genes = set()
for tissue in zvals.columns:
    tvals = zvals[[tissue]].copy()
    tvals["gene"] = [item.split("_")[0] for item in tvals.index]
    tvals["ID"] = tvals.index
    tdata = tvals.groupby("gene", as_index=False).agg({tissue: max})
    tdata.columns = ["gene","max"]
    tvals = pd.merge(tvals, tdata, on=["gene"])
    tvals = tvals[tvals[tissue]==tvals["max"]]
    tvals_sig = tvals[tvals["max"]>=ZTHRESH]
    print("%s: %s gene-level eSTRs"%(tissue, tvals_sig.shape[0]))
    all_estrs = all_estrs.union(set(tvals_sig["ID"]))
    all_genes = all_genes.union(set(tvals_sig["gene"]))
    tvals_sig.index = tvals_sig["ID"]
    tvals_sig[[tissue]].to_csv(os.path.join(workdir, "output-%s"%prefix, "sig-bytissue", "%s-estrs.tsv"%tissue), sep="\t", index=True)

# Summarize - overall
print("Sig genes (global best per gene): %s"%sig.shape[0])
print("Sig genes (best per tissue for each gene): %s eSTRs, %s unique genes"%(len(all_estrs), len(all_genes)))
print("Sig units: %s"%(zmax[zmax["zmax"]>=ZTHRESH].shape[0]))

