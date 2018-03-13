#!/usr/bin/env python
"""
Usage:
./fisher_exact.py <table>
table is tab delim with:
category, num_AA, num_AB, num_BA, numBB
Output table with p-value column
"""

import sys
import pandas as pd
import scipy.stats

try:
    tfile = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

table = pd.read_csv(tfile, sep="\t", names=["category", "AA","AB","BA","BB"])

def GetFisher(x):
    t = [[x["AA"], x["AB"]],[x["BA"],x["BB"]]]
    x = scipy.stats.fisher_exact(t)
    return ",".join(map(str, x))

res = table.apply(GetFisher, 1)
table["oddsratio"] = res.apply(lambda x: x.split(",")[0])
table["pval"] = res.apply(lambda x: x.split(",")[1])

table.to_csv(sys.stdout, sep="\t", index=False)
