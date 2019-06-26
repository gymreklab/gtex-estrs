#!/usr/bin/env python3

import argparse
import glob
import numpy as np
import pandas as pd
import sys

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
    anova.ix[anova["SNP"].apply(str) == "nan","anova.pval"] = -1*np.inf # no SNP available
    def GetSNP(x):
        try:
            return int(x.split(":")[1])
        except: return -1
    anova["mashr.top.snp"] = anova["SNP"].apply(GetSNP)
    anova["anova.best.snp"] = anova["SNP"]
    return anova[["gene","str.start","mashr.top.snp","anova.pval","anova.best.snp"]]

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
    caviar = pd.concat(dfl, axis=0, ignore_index=True) #gene    top_snp top_snp_score   str     str.score       str.rank        num.snps
    caviar["str.start"] = caviar["str"].apply(lambda x: int(x.split(":")[-1]))
    caviar["caviar.str.score"] = caviar["str.score"]
    caviar["caviar.str.rank"] = caviar["str.rank"]
    caviar["caviar.topsnp"] = caviar["top_snp"]
    caviar["caviar.topsnp.score"] = caviar["top_snp_score"]
    caviar["caviar.nsnps"] = caviar["num.snps"]
    return caviar[["gene","str.start","caviar.str.score","caviar.str.rank","caviar.topsnp","caviar.topsnp.score","caviar.nsnps"]]

def LoadGeneAnnot(geneannotfile):
    annot = pd.read_csv(geneannotfile, sep="\t")
    return annot[["gene","gene.name","gene.strand"]]

def CheckCols(data, cols):
    for col in cols:
        numNA = sum(data[col].apply(str)=="nan")
        if numNA > 0:
            sys.stderr.write(str(data[data[col].apply(str)=="nan"]))
            sys.stderr.write("\nError: column %s has %s NA values\n"%(col, numNA))
            return False
    return True

def CheckTable(data):
    # Check no NAs in universal columns
    if not CheckCols(data, ["gene","chrom","str.start", \
                            "mashr.beta", \
                            "str.motif.forward","str.motif.reverse"]): return False
    # If we have caviar, we should have anova and vice versa
    caviarNA = np.isnan(data["caviar.str.score"])
    anovaNA = (data["anova.best.snp"]=="NA")
    if not all([caviarNA[i] == anovaNA[i] for i in range(data.shape[0])]):
        sys.stderr.write(str(data[anovaNA != caviarNA][["chrom","str.start","gene","caviar.str.score","anova.pval"]]))
        sys.stderr.write("\nERROR: Caviar and anova don't match\n")
        return False
    # If we have mashr.significant, we should have anova and caviar (unless they exploded?)
    mashrSig = data["mashr.significant"]
    if not all ([mashrSig[i] != caviarNA[i] for i in range(data.shape[0])]):
        sys.stderr.write(str(data[mashrSig[i]==caviarNA[i]]))
        sys.stderr.write("\nERROR: Missing CAVIAR for mashr.significant loci")
        return False
    if not all ([mashrSig[i] != anovaNA[i] for i in range(data.shape[0])]):
        sys.stderr.write(str(data[mashrSig[i]==anovaNA[i]]))
        sys.stderr.write("\nERROR: Missing anova for mashr.significant loci")
        return False
    return True

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
    if not CheckCols(data, ["chrom","gene","str.start"]): sys.exit(1)

    # Load linreg
    linreg = LoadLinreg(args.linreg)
    data = pd.merge(data, linreg, on=["gene","str.start"], how="outer")
    if not CheckCols(data, ["chrom","gene","str.start"]): sys.exit(1)
    
    # Load ANOVA
    anova = LoadAnova(args.anova)
    data = pd.merge(data, anova, on=["gene", "str.start"], how="outer")
    if not CheckCols(data, ["chrom","gene","str.start"]): sys.exit(1)

    # Load CAVIAR
    caviar = LoadCaviar(args.caviar)
    data = pd.merge(data, caviar, on=["gene","str.start"], how="outer")
    if not CheckCols(data, ["chrom","gene","str.start"]): sys.exit(1)

    # Load gene annotations
    annot = LoadGeneAnnot(args.geneannot)
    data = pd.merge(data, annot, on=["gene"])

    # Check table and output
    if CheckTable(data):
        data.sort_values("caviar.str.score", ascending=False).to_csv(args.out, sep="\t", index=False)
