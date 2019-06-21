#!/usr/bin/env python3

"""
gene
gene.name
gene.strand
chrom
str.start (hipref)
str.end (hipref)
str.motif.forward (hipref)
str.motif.reverse (hipref)
str.id
linreg.beta (linreg)
linreg.beta.se (linreg)
linreg.p.wald (linreg)
linreg.top.str (linreg)
linreg.n.miss (linreg)
mashr.beta (mashr)
mashr.beta.se (mashr)
mashr.top.str (mashr)
mashr.top.snp (mashr - get from ANOVA)
anova.pval (anova)
anova.qval (anova)
caviar.str.score (caviar)
caviar.topsnp (caviar)
caviar.topsnp.score (caviar)
caviar.str.rank (caviar)
significant (mashr - top.str)
"""

import argparse
import pandas as pd

def LoadMashr(mashr_beta_file, mashr_se_file):
    mashr_betas = pd.read_csv(args.mashr_beta, sep="\t", index_col=0)
    mashr_ses = pd.read_csv(args.mashr_se, sep="\t", index_col=0)
    mashr_genes = [item.split("_")[0] for item in mashr_betas.index]
    mashr_chroms = [item.split("_")[1] for item in mashr_betas.index]
    mashr_pos = [int(item.split("_")[2]) for item in mashr_betas.index]
    mashr = pd.DataFrame({"mashr.beta": mashr_betas[args.tissue], \
                          "mashr.beta.se": mashr_ses[args.tissue], \
                          "chrom": mashr_chroms, \
                          "gene": mashr_genes, \
                          "str.start": mashr_pos})
    mashr["Z"] = abs(mashr["mashr.beta"]/mashr["mashr.beta.se"])
    mashr_best = mashr.groupby("gene", as_index=False).agg({"Z": max})
    mashr_best.columns = ["gene","maxZ"]
    mashr = pd.merge(mashr, mashr_best, on=["gene"])
    mashr["mashr.top.str"] = (mashr["Z"]==mashr["maxZ"])
    return mashr[["chrom","gene","str.start","mashr.beta","mashr.beta.se","mashr.top.str"]]

def LoadLinreg(linregfile):
    linreg = pd.read_csv(linregfile, sep="\t")
    linreg = linreg[["gene","str.start","beta","beta.se","p.wald","n.miss"]]
    linreg.columns = ["gene","str.start","linreg.beta","linreg.beta.se","linreg.pval","linreg.n.miss"]
    return linreg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CAVIAR on GTEx data")
    parser.add_argument("--hipref", help="HipSTR reference file with stranded motif info", type=str, required=True)
    parser.add_argument("--tissue", help="Tissue to use", type=str, required=True)
    parser.add_argument("--linreg", help="STR linear regression results file", type=str, required=True)
    parser.add_argument("--mashr-beta", help="File with mashR posterior betas", type=str, required=True)
    parser.add_argument("--mashr-se", help="File with mashR posterior std errs", type=str, required=True)
    args = parser.parse_args()

    # Load HipSTR
    hipref = pd.read_csv(args.hipref, sep="\t", names=["chrom","str.start","str.end","period","str.motif.forward","str.motif.reverse"])

    # Load mashR
    mashr = LoadMashr(args.mashr_beta, args.mashr_se)
    data = pd.merge(mashr, hipref[["chrom","str.start","str.end","str.motif.forward","str.motif.reverse"]], on=["chrom","str.start"])

    # Load linreg
    linreg = LoadLinreg(args.linreg)
    data = pd.merge(data, linreg, on=["gene","str.start"], how="outer")
    print(data.head())
