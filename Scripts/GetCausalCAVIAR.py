import argparse
import math
import numpy as np
import os
import pandas as pd
import random
import shutil
import sys
import gzip

"""
Quantify the probability of a variant to be causal while allowing with arbitrary number of causal variants
"""

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())
        
def z(vals):
    vals['z.score']=vals['beta']/vals['beta.se']
    return vals['z.score']

def WriteCorrTable(indexed_genotypes):
    """ generate correlation table using normalized genotype
      _1_ ... _n_
    1| 1  ... C1n
    .|    ...
    n|Cn1 ... Cnn=1
    """
    G=indexed_genotypes.transpose()
    variants = list(G.columns)
    frames=[]
    len(set(variants))
    for V1 in variants:
        COV=[]
        for V2 in variants:
            X=G[V1].replace('None', np.nan).astype(float)
            Y=G[V2].replace('None', np.nan).astype(float)
            COV.append(X.corr(Y))
        frames.append(COV)
    return pd.DataFrame(frames,columns=variants, index=variants) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Linreg/LMM simulations with real data")
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    parser.add_argument("--exprannot", help="Expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict analysis to this chromosome", type=str, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--snpgt", help="File with normalized SNP genotypes", type=str, required=True)
    parser.add_argument("--out", help="Write results to this file", type=str, required=True)
    parser.add_argument("--tmpdir", help="Use this directory for temporary files", type=str, default="/tmp")
    parser.add_argument("--linreg_snp", help="File with snp linear regression output", type=str, required=True)
    parser.add_argument("--linreg_str", help="File with str linear regression output", type=str, required=True)
    parser.add_argument("--debug", help="Print debug status messages", action="store_true")
    parser.add_argument("--restrict_to_estrs", help="Restrict the analysis to eSTRs genes", type=str, required=False)
    parser.add_argument("--distfromgene", help="Look at STRs/SNPs within this distance of gene boundaries", type=int, required=True)

    args = parser.parse_args()
    EXPRFILE = args.expr
    EXPRANNOTFILE = args.exprannot
    CHROM = args.chrom
    if "chr" not in str(CHROM): CHROM="chr%s"%CHROM
    DISTFROMGENE = args.distfromgene
    STRGTFILE = args.strgt
    SNPGTFILE = args.snpgt
    OUTFILE = args.out
    ESTRGENESFILE = args.restrict_to_estrs
    REG_SNPs = args.linreg_snp
    REG_STRs = args.linreg_str
    TMPDIR = args.tmpdir
    DEBUG = args.debug

    # Load expression
    PROGRESS("\nLoad expression", printit=DEBUG)
    expr = pd.read_csv(EXPRFILE)

    # Load annotation
    PROGRESS("Load annotation", printit=DEBUG)
    expr_annot = pd.read_csv(EXPRANNOTFILE)
    expr_annot.index = expr_annot["probe.id"].values
    expr_annot = expr_annot.loc[list(expr.columns)].dropna() 
    expr_annot = expr_annot[expr_annot["gene.chr"] == CHROM]

    # Load strs Regression
    PROGRESS("\nLoad strs regression", printit=DEBUG)
    strs = pd.read_csv(REG_STRs, sep="\t")
    strs = strs.loc[strs['chrom']==CHROM]

    # Load snps regression
    PROGRESS("\nLoad snps regression", printit=DEBUG)
    snps = pd.read_csv(REG_SNPs, sep="\t")
    snps = snps.loc[snps['chrom']==CHROM]

    #Load SNP genotypes
    PROGRESS("Load SNPs", printit=DEBUG)
    snpgt = pd.read_csv(SNPGTFILE, sep="\t")
    snpgt = snpgt.loc[snpgt['chrom']==CHROM]

    # Load STR genotypes
    PROGRESS("Load STRs", printit=DEBUG)
    strgt = pd.read_csv(STRGTFILE, sep="\t")
    strgt = strgt.loc[strgt['chrom']==CHROM]

    # Restrict to STR samples
    PROGRESS("Restrict to STRs samples", printit=DEBUG)
    str_samples = list(set(strgt.columns[2:].values).intersection(set(snpgt.columns[2:].values)))
    expr = expr.loc[str_samples,:]
    snpgt = snpgt[["chrom","start"] + str_samples]
    snpgt.index = list(snpgt["start"].apply(lambda x: "SNP_%s"%int(x)))
    strgt = strgt[["chrom","start"] + str_samples]
    strgt.index = list(strgt["start"].apply(lambda x: "STR_%s"%int(x)))
    samples_to_keep = str_samples

    # Load eSTR results
    PROGRESS("Restrict to eSTR genes only", printit=DEBUG)
    if ESTRGENESFILE is not None:
        estr_genes = pd.read_csv(ESTRGENESFILE, sep="\t")
        Genes = estr_genes.loc[estr_genes['qvalue']<=0.1]['gene']  # estrs at 10%FDR
        expr_annot = expr_annot.loc[expr_annot['gene.id'].isin(list(Genes))]

# For each gene, get all cis-variants and the best STR
    for i in range(expr_annot.shape[0]):
        gene=expr_annot.index.values[i]
        ensgene = expr_annot["gene.id"].values[i]
        genedir=TMPDIR+"/%s"%gene
        if not os.path.exists(genedir):
            os.mkdir(genedir)
        PROGRESS("Getting data for %s"%gene, printit=DEBUG)
        start = expr_annot["gene.start"].values[i]
        end = expr_annot["gene.stop"].values[i]
    # Pull out cis SNPs
        PROGRESS("Getting cis SNPs for %s"%gene)
        cis_snps = snps[(snps["str.start"] >= (start-DISTFROMGENE)) & (snps["str.start"] <= (end+DISTFROMGENE))]
        cis_snps = cis_snps.loc[cis_snps['gene']==ensgene]
        cis_snps.index = cis_snps["str.start"].apply(lambda x: "SNP_%s"%int(x))
        L=list(cis_snps.index)
    # Pull out most significant STR
        PROGRESS("Getting most significant cis STR for %s"%gene)
        best_str_start = strs[strs["gene"]==ensgene].sort("p.wald")["str.start"].values[0]
        cis_strs = strs.loc[strs['gene']==ensgene]
        cis_strs.index = list(cis_strs['str.id'])
        try:
            del cis_strs['ID']
        except:
            pass
        cis_snps.loc['STR_'+str(best_str_start)] = list(cis_strs.loc['STR_'+str(best_str_start)])
    # Make z file data
        Z = z(cis_snps[['beta','beta.se']])
        Z.to_csv(genedir+'/ZFILE', sep='\t',header=None)
    # Make LD file
        genotypes = snpgt.loc[L]
        genotypes.loc['STR_'+str(best_str_start)] = list(strgt.loc['STR_'+str(best_str_start)])
        del genotypes['chrom']
        del genotypes['start']
        Matrix = WriteCorrTable(genotypes)
        Matrix.to_csv(genedir+'/LDFILE', sep='\t',header=None)
        PROGRESS("Matrix of corr was sent to file for %s"%gene)
    #Run caviar
        caviar_cmd = "/usr/local/bin/CAVIAR -l %s -z %s -o %s/caviar -r -c 1 -f 1 "%(genedir +"/LDFILE", genedir+"/ZFILE", genedir)
        os.system(caviar_cmd)
    # Output results
    
    # Need to add code for merging results into one file. At this point caviar output for each gene
    # resides as "caviar" in the dir of gene name within tmp director. Will need to merge them
    





