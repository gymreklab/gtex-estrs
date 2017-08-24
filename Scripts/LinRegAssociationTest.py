import argparse
import math
import numpy as np
import os
import pandas as pd
import random
import shutil
import statsmodels.api as sm
import sys

EXPRFILE = None
EXPRANNOTFILE = None
DISTFROMGENE = None
STRGTFILE = ""
OUTFILE = ""
SCRIPTDIR = "/san/melissa/workspace/str-qtl/lmm"
TMPDIR = "/tmp/"
CHECKCHROM = False
PERMUTE_EXPR = False
DEBUG = False
MINSAMPLES = 0

def PROGRESS(msg):
    sys.stderr.write("%s\n"%msg.strip())

def ZNorm(vals):
    m = np.mean(vals)
    sd = math.sqrt(np.var(vals))
#    print m, '  ', sd
    if sd == 0: return None
    return [(item-m)/sd for item in vals]

def LinearRegression(X, Y, norm=False, minsamples=0):
    """
    Perform linear regression, return beta, beta_se, p
    """
#    print 'This is X \t', X, type(X)
#    print 'And this Y\t', Y, type(Y)
    if norm:
        X = ZNorm(X)
        Y = ZNorm(Y)
        if X is None or Y is None: return None, None, None
        if np.var(X)==0: return None, None, None
        if len(X) <= minsamples: return None, None, None
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X)
    res_ols = mod_ols.fit()
    pval = res_ols.pvalues[1]
    print 'P-value: ', pval 
    slope = res_ols.params[1]
    err = res_ols.bse[1]
#    print 'slope: ', slope, ' Err: ', err
    return slope, err, pval

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get datasets for LMM assocation tests with STRs and SNPs")
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    parser.add_argument("--exprannot", help="Expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict analysis to this chromosome", type=str, required=True)
    parser.add_argument("--distfromgene", help="Look at STRs/SNPs within this distance of gene boundaries", type=int, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--out", help="Write data files to this file", type=str, required=True)
    parser.add_argument("--scriptdir", help="Directory containing scripts", type=str, required=False)
    parser.add_argument("--checkchrom", help="Only load exons for relevant chromosome", action="store_true")
    parser.add_argument("--tmpdir", help="Tmp directory", type=str)
    parser.add_argument("--permute", help="Permute expression values", action="store_true")
    parser.add_argument("--norm", help="Normalize genotypes before doing association", action="store_true")
    parser.add_argument("--min-samples", help="Require data for this many samples", type=int, default=0)
    parser.add_argument("--debug", help="Print debug info", action="store_true")
    args = parser.parse_args()
    EXPRFILE = args.expr
    EXPRANNOTFILE = args.exprannot
    CHROM = args.chrom
    if "chr" not in str(CHROM): CHROM="chr%s"%CHROM
    DISTFROMGENE = args.distfromgene
    STRGTFILE = args.strgt
    OUTFILE = args.out
    NORM = args.norm
    MINSAMPLES = args.min_samples
    if args.scriptdir is not None: SCRIPTDIR = args.scriptdir
    if args.checkchrom: CHECKCHROM = True
    if args.tmpdir is not None: TMPDIR = args.tmpdir
    if args.permute: PERMUTE_EXPR = True
    if args.debug: DEBUG = True

    # Load expression values
    PROGRESS("Load expression")
    if CHECKCHROM:
        x = list(pd.read_csv(EXPRFILE, nrows=1).columns.values)
        x = [item for item in x if item == "Unnamed: 0" or CHROM in item]
        expr = pd.read_csv(EXPRFILE, usecols=x) # reading in all the columns takes way too much memory,
                                                # pull out only ones we need
    else:
        expr = pd.read_csv(EXPRFILE)
    if "Unnamed: 0" in expr.columns:
	expr.index = expr["Unnamed: 0"].values
    	expr = expr.drop("Unnamed: 0", 1)
#    print expr

    # Load expression annotation
    PROGRESS("Load annotation")
    expr_annot = pd.read_csv(EXPRANNOTFILE)
#    print expr_annot.head(),'1'
    expr_annot.index = expr_annot["probe.id"].values
#    print expr_annot.head(), '2'
#    print expr_annot.columns, '2'
    expr_annot = expr_annot.loc[[item for item in expr.columns if item in expr_annot.index],:]
#    print expr_annot.head(), '3'
    expr_annot = expr_annot[expr_annot["gene.chr"] == CHROM]
#    print expr_annot.head(),'4'
    # Load STR genotypes
    PROGRESS("Load STRs")
    strgt = pd.read_csv(STRGTFILE, sep="\t")
    print strgt.head()
    strgt = strgt[strgt["chrom"] == CHROM]   ; print "Restrict to chrom"
    # Restrict to STR samples
    str_samples = list(set(strgt.columns[2:].values))
#    print expr_annot.head()

#    print len(str_samples), expr.shape
    samples_to_remove = []
    for item in str_samples:
        if item not in expr.index: samples_to_remove.append(item) #str_samples.remove(item)
    for item in samples_to_remove: str_samples.remove(item)
    expr = expr.loc[str_samples,:]
#    print len(str_samples), '   ', expr.shape

#    print expr
    f = open(OUTFILE, "w")
    f.write("\t".join(["gene","chrom","str.id","str.start","n.miss","allele1.dummy","allele2.dummy","af.dummy","beta","beta.se","lambda.remel","p.wald"])+"\n")
    # For each gene:
    # Pull out STRs within distance of gene ends
    # For each STRxgene pair, get expr, str for samples with data and do linreg
    for i in range(expr_annot.shape[0]):
        gene = expr_annot.index.values[i]
        PROGRESS("Getting data for %s"%gene)
        start = expr_annot["gene.start"].values[i]
        end = expr_annot["gene.stop"].values[i]
        cis_strs = strgt[(strgt["start"] >= (start-DISTFROMGENE)) & (strgt["start"] <= (end+DISTFROMGENE))]
        if PERMUTE_EXPR:
            expr[gene] = random.sample(list(expr[gene].values), expr.shape[0])
        y = pd.DataFrame({"expr":list(expr.loc[:, gene])})
        y.index = str_samples
#	print cis_strs
#	print y
        for j in range(cis_strs.shape[0]):
            # Get STR data
            locus_str = cis_strs.iloc[[j],:][str_samples].transpose()
            locus_str.index = str_samples
            locus_str.columns = ["STR_%s"%(cis_strs["start"].values[j])]
            samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if str(locus_str.iloc[:,0].values[k]) != "None"]
            # Get subsets
            locus_str = locus_str.loc[samples_to_keep,:]
            locus_y = y.loc[samples_to_keep,:]
            print locus_str
            # Run regression
#           print ' Regression', gene
            beta, beta_se, p = LinearRegression(map(float,locus_str.iloc[:,0].values), locus_y["expr"].values, norm=NORM, minsamples=MINSAMPLES)
#	    print 'Pvalue: ', p
            # Write output
            if beta is not None:
                str_start = cis_strs["start"].values[j]
                f.write("\t".join(map(str, [gene, CHROM, "STR_%s"%str_start, str_start, len(str_samples)-locus_str.shape[0], "A", "G", 0, beta, beta_se, -1, p]))+"\n")
    f.close()
