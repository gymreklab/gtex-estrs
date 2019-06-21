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
import glob
import numpy as np
import pandas as pd

ZTHRESH = 3

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
    mashr["mashr.significant"] = mashr.apply(lambda x: x["mashr.top.str"] and abs(x["Z"])>=ZTHRESH, 1)
    mashr["str.start"] = mashr["str.start"].apply(int)
    return mashr[["chrom","gene","str.start","mashr.beta","mashr.beta.se","mashr.top.str","mashr.significant"]]

def LoadLinreg(linregfile):
    linreg = pd.read_csv(linregfile, sep="\t")
    linreg = linreg[["gene","str.start","beta","beta.se","p.wald","n.miss"]]
    linreg.columns = ["gene","str.start","linreg.beta","linreg.beta.se","linreg.pval","linreg.n.miss"]
    linreg["str.start"] = linreg["str.start"].apply(int)
    return linreg

def LoadAnova(anovafiles):
    allfiles = glob.glob(anovafiles)
    dfl = []
    for fname in allfiles:
        df = pd.read_csv(fname, sep="\t")
        df.columns = ["gene","STR","SNP","anova.pval"]
        dfl.append(df)
    anova = pd.concat(dfl, axis=0, ignore_index=True)
    anova["str.start"] = anova["STR"].apply(lambda x: int(x.split(":")[1]))
    def GetSNP(x):
        try:
            return int(x.split(":")[1])
        except: return -1
    anova["mashr.top.snp"] = anova["SNP"].apply(GetSNP)
    return anova[["gene","str.start","mashr.top.snp","anova.pval"]]

def LoadHipref(hipreffile):
    hipref = pd.read_csv(hipreffile, sep="\t", names=["chrom","str.start","str.end","period","str.motif.forward","str.motif.reverse"])
    hipref["str.start"] = hipref["str.start"].apply(int)
    hipref["str.end"] = hipref["str.end"].apply(int)
    return hipref

def LoadCaviar(caviarfiles):
    allfiles = glob.glob(caviarfiles)
    dfl = []
    for fname in allfiles:
        df = pd.read_csv(fname, sep="\t")
        dfl.append(df)
    caviar = pd.concat(dfl, axis=0, ignore_index=True)
    caviar["str.start"] = caviar["top_str"].apply(lambda x: int(x.split(":")[-1]))
    caviar["caviar.str.score"] = caviar["top.str.score"]
    caviar["caviar.str.rank"] = caviar["str.rank"]
    caviar["caviar.topsnp"] = caviar["top_snp"]
    caviar["caviar.topsnp.score"] = caviar["top_snp_score"]
    caviar["caviar.nsnps"] = caviar["num.snps"]
    return caviar[["gene","str.start","caviar.str.score","caviar.str.rank","caviar.topsnp","caviar.topsnp.score","caviar.nsnps"]]

def LoadGeneAnnot(geneannotfile):
    annot = pd.read_csv(geneannotfile, sep="\t")
    return annot[["gene","gene.name","gene.strand"]]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CAVIAR on GTEx data")
    parser.add_argument("--hipref", help="HipSTR reference file with stranded motif info", type=str, required=True)
    parser.add_argument("--tissue", help="Tissue to use", type=str, required=True)
    parser.add_argument("--linreg", help="STR linear regression results file", type=str, required=True)
    parser.add_argument("--mashr-beta", help="File with mashR posterior betas", type=str, required=True)
    parser.add_argument("--mashr-se", help="File with mashR posterior std errs", type=str, required=True)
    parser.add_argument("--out", help="Output file name", type=str, required=True)
    parser.add_argument("--anova", help="Path to ANOVA files", type=str, required=True)
    parser.add_argument("--caviar", help="Path to CAVIAR files", type=str, required=True)
    parser.add_argument("--geneannot", help="Path to gene annotations", type=str, required=True)
    args = parser.parse_args()

    # Load HipSTR
    hipref = LoadHipref(args.hipref)

    # Load mashR
    mashr = LoadMashr(args.mashr_beta, args.mashr_se)
    data = pd.merge(mashr, hipref[["chrom","str.start","str.end","str.motif.forward","str.motif.reverse"]], on=["chrom","str.start"])

    # Load linreg
    linreg = LoadLinreg(args.linreg)
    data = pd.merge(data, linreg, on=["gene","str.start"], how="outer")
    
    # Load ANOVA
    anova = LoadAnova(args.anova)
    data = pd.merge(data, anova, on=["gene", "str.start"], how="outer")

    # Load CAVIAR
    caviar = LoadCaviar(args.caviar)
    data = pd.merge(data, caviar, on=["gene","str.start"], how="outer")

    annot = LoadGeneAnnot(args.geneannot)
    data = pd.merge(data, annot, on=["gene"])

    # Output
    data.sort_values("caviar.str.score", ascending=False).to_csv(args.out, sep="\t", index=False)
