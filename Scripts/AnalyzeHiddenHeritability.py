#!/usr/local/bin/python2.7

import argparse
import math
import numpy as np
import os
import pandas as pd
import sys
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

"""
Analyze h2 gain in adding eSTRs vs. eSNPs
"""

def ZNorm(vals):
    vals = list(vals)
    #print vals
    vals = [float(x) for x in vals]
    m = np.mean(vals)
    sd = math.sqrt(np.var(vals))
    print m,'=============', sd
    if sd == 0: return None
    return [(item-m)*1.0/sd for item in vals]

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())
DISTFROMGENE = 100000

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Linreg/LMM simulations with real data")
    parser.add_argument("--estrs", help="Estr results", type=str, required=True)
    parser.add_argument("--esnps", help="Esnp results", type=str, required=True)
    parser.add_argument("--str-fdr-method", help="FDR method (STRs)", type=str, required=False, default="qval.gene")
    parser.add_argument("--str-fdr-threshold", help="Cutoff for FDR (STRs)", type=float, required=False, default=1)
    parser.add_argument("--snp-fdr-method", help="FDR method (SNPs)", type=str, required=False, default="qval.gene")
    parser.add_argument("--snp-fdr-threshold", help="Cutoff for FDR (SNPs)", type=float, required=False, default=1)
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    parser.add_argument("--exprannot", help="Expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict analysis to this chromosome", type=str, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--snpgt", help="File with normalized SNP genotypes", type=str, required=True)
    parser.add_argument("--nohead", help="Don't print header", action="store_true")
    parser.add_argument("--debug", help="Print debug status messages", action="store_true")
    args = parser.parse_args()
    ESTR_RESULTS_FILE = args.estrs
    ESNP_RESULTS_FILE = args.esnps
    SNP_FDR_METHOD = args.snp_fdr_method
    SNP_FDR_THRESHOLD = args.snp_fdr_threshold
    STR_FDR_METHOD = args.str_fdr_method
    STR_FDR_THRESHOLD = args.str_fdr_threshold
    EXPRFILE = args.expr
    EXPRANNOTFILE = args.exprannot
    CHROM = args.chrom
    if "chr" not in str(CHROM): CHROM="chr%s"%CHROM
    STRGTFILE = args.strgt
    SNPGTFILE = args.snpgt
    NOHEAD = args.nohead
    DEBUG = args.debug

    # Load expression and annotation
    PROGRESS("Load expression", printit=DEBUG)
    expr = pd.read_csv(EXPRFILE)
    print expr.shape
#    expr.index = expr["Unnamed: 0"].values
#    expr = expr.drop("Unnamed: 0", 1)
    PROGRESS("Load annotation", printit=DEBUG)
    expr_annot = pd.read_csv(EXPRANNOTFILE)
    expr_annot.index = expr_annot["probe.id"].values
    expr_annot = expr_annot.loc[expr.columns[map(lambda x: x in expr_annot.index, expr.columns)],:]
    expr_annot = expr_annot[expr_annot["gene.chr"] == CHROM]
    print expr_annot.shape
    
    # Load SNP genotypes
    PROGRESS("Load SNPs", printit=DEBUG)
    snpgt = pd.read_csv(SNPGTFILE, sep="\t")
    print snpgt.shape

    # Load STR genotypes
    PROGRESS("Load STRs", printit=DEBUG)
    strgt = pd.read_csv(STRGTFILE, sep="\t")
    print strgt.shape

    # Restrict to STR samples
    str_samples = list(set(strgt.columns[2:].values).intersection(set(snpgt.columns[2:].values)))
    expr = expr.loc[str_samples,:]
    snpgt = snpgt[["chrom","start"] + str_samples]
    strgt = strgt[["chrom","start"] + str_samples]

    # Load eQTL results
    PROGRESS("Load eQTL results", printit=DEBUG)
    estr_results = pd.read_csv(ESTR_RESULTS_FILE, sep="\t")
    esnp_results = pd.read_csv(ESNP_RESULTS_FILE, sep="\t")
    print 'eqtl outputs loaded ', estr_results.shape, '\t',esnp_results.shape

    # Print output header
    if not NOHEAD:
        print ",".join(["chrom","gene","str.start","numsnps","numsamples","r2_str","r2_snp","r2_snpstr","anova_pval","estr_fdr","esnp_fdr","delta_bic","delta_aic"])

    # For each gene in results:
    # Pull out passing eSNPs and eSTRs
    # Build STR, SNP, and SNPSTR model
    # Get r2 for each model and run Anova + BIC
    print CHROM,'---fdr-method----', STR_FDR_METHOD, '---treshold---',STR_FDR_THRESHOLD 
    
    estr_results = estr_results[(estr_results["chrom"]==CHROM) & (estr_results[STR_FDR_METHOD]<=STR_FDR_THRESHOLD)]
    print estr_results.columns
    
    if STR_FDR_METHOD == "qval.gene":
        estr_results = estr_results[estr_results["best_str"]==1]
        print estr_results.columns,'*******\n'
        
        
    esnp_results = esnp_results[(esnp_results["chrom"]==CHROM) & (esnp_results[SNP_FDR_METHOD]<=SNP_FDR_THRESHOLD)]
    print esnp_results.shape , '------\n', esnp_results.columns
    
    
    if SNP_FDR_METHOD == "qval.gene":
        esnp_results = esnp_results[esnp_results["best_str"]==1]
    
    allgenes = set(estr_results['gene'])          #set(estr_results["gene"].values)
    for ensgene in allgenes:
        estr_fdr = min(estr_results[estr_results["gene"]==ensgene][STR_FDR_METHOD])
        esnp_fdr = min(esnp_results[esnp_results["gene"]==ensgene][SNP_FDR_METHOD])
        probes = expr_annot[expr_annot["gene.id"]==ensgene]["probe.id"].values
        if len(probes) != 1:
            PROGRESS("ERROR: no unique probe for for %s"%ensgene)
            continue
        else: probe = probes[0]
        PROGRESS("Getting data for %s (%s)"%(ensgene, probe), printit=DEBUG)
        # Get eSTRs
        estrs = estr_results[(estr_results["gene"]==ensgene)]
        # get eSNPs
        esnps = esnp_results[(esnp_results["gene"]==ensgene)]
        # Get genotypes and expression values
        y = pd.DataFrame({"expr":list(expr.loc[:,probe])})
        y.index = str_samples
        cis_snps = snpgt[snpgt["start"].apply(lambda x: x in list(esnps["str.start"].values))]
        cis_strs = strgt[strgt["start"].apply(lambda x: x in list(estrs["str.start"].values))]
        # Run one analysis per STR
        for i in range(cis_strs.shape[0]):
            # Get STR data
            locus_str = cis_strs.iloc[[i],:][str_samples].transpose()
            locus_str.index = str_samples
            locus_str.columns = ["STR_%s"%(cis_strs["start"].values[i])]
            samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if str(locus_str.iloc[:,0].values[k]) != "None"]
            # Get subsets
            locus_str = locus_str.loc[samples_to_keep,:]
            locus_snp = cis_snps[str_samples].transpose()
            locus_snp = locus_snp.loc[samples_to_keep,:]
            locus_y = y.loc[samples_to_keep,:]
            print ZNorm(locus_y["expr"].values), '\n znorm or expr itself?\n'
            print locus_y["expr"].values
            # Get data frame with relevant variables
            d = {"STR": ZNorm(locus_str.iloc[:,0].apply(float).values), "expr": ZNorm(locus_y["expr"].values)}
            for k in range(locus_snp.shape[1]):
                snps = ZNorm(locus_snp.iloc[:,k].values)
                print snps
                if snps is not None:
                    d["SNP%s"%k] = snps
            genedata = pd.DataFrame(d)
            # Run regression with SNPs
            if locus_snp.shape[1] > 0 and "SNP0" in genedata.columns:
                snpstring = "+".join([item for item in list(genedata.columns) if "SNP" in item])
                formula_snpstr = "expr ~ STR + " + snpstring
                formula_snp = "expr ~ " + snpstring
                formula_str = "expr ~ STR"
                print genedata
                lm_snpstr = ols(formula_snpstr, genedata).fit()
                lm_snp = ols(formula_snp, genedata).fit()
                lm_str = ols(formula_str, genedata).fit()
                bic_snpstr = lm_snpstr.bic
                bic_snp = lm_snp.bic
                aic_snpstr = lm_snpstr.aic
                aic_snp = lm_snp.aic
                snpstr_rsq = lm_snpstr.rsquared
                snp_rsq = lm_snp.rsquared
                str_rsq = lm_str.rsquared
                anova_results = anova_lm(lm_snp, lm_snpstr)
                pval = anova_results["Pr(>F)"].values[1]
                print ",".join(map(str, [CHROM, ensgene, cis_strs["start"].values[i], locus_snp.shape[1],\
                                              locus_y.shape[0], str_rsq, snp_rsq, snpstr_rsq, pval, estr_fdr, esnp_fdr,\
                                             bic_snp-bic_snpstr, aic_snp-aic_snpstr]))
            else:
                formula_str = "expr ~ STR"
                lm_str = ols(formula_str, genedata).fit()
                str_rsq = lm_str.rsquared
                print ",".join(map(str, [CHROM, ensgene, cis_strs["start"].values[i], locus_snp.shape[1],\
                                              locus_y.shape[0], str_rsq, 0, str_rsq, -1, estr_fdr, -1,\
                                             0, 0]))
                