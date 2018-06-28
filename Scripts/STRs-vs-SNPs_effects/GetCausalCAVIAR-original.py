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
Quantify the probability of a variant to be causal. We may later allow for arbitrary number of causal variants
"""

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())
        
def MakeZScoreTable(vals):
    try:
        Zscore=vals['beta']/vals['beta.se']
    except:
        Zscore=None
    return Zscore

def WriteCorrTable(indexed_genotypes):
    """ generate correlation table using normalized genotype
      _1_ ... _n_
    1| 1  ... C1n
    .|    ...
    n|Cn1 ... Cnn=1
    """
    G=indexed_genotypes.transpose()
    variants = list(G.columns)
    CMat=[]
    print '\t\t\t**', len(variants)
    for V1 in variants:
        COV=[]
        for V2 in variants:
            X=G[V1].replace('None', np.nan).astype(float)
            Y=G[V2].replace('None', np.nan).astype(float)
            #COV.append(X.corr(Y))
	    if X.corr(Y) is np.nan:    #### 
		print X, list(X),'\n**',Y, list(Y)
		COV.append(0.0) #For missing LD we assume non linear corr (undetermined LD)
	    else:
		COV.append(X.corr(Y))
        CMat.append(COV)
    return pd.DataFrame(CMat,columns=variants, index=variants) 

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
    parser.add_argument("--restrict_to_estrs", help="Restrict the analysis to eSTRs genes from this file", type=str, required=False)
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
    if not os.path.exists(TMPDIR):
        os.mkdir(TMPDIR)

    # Load expression
    PROGRESS("\nLoad expression", printit=DEBUG)
    expr = pd.read_csv(EXPRFILE)

    # Load annotation
    PROGRESS("Load annotation", printit=DEBUG)
    expr_annot = pd.read_csv(EXPRANNOTFILE)
    expr_annot.index = expr_annot["probe.id"].values
    expr_annot = expr_annot.reindex(list(expr.columns))
    expr_annot = expr_annot.dropna() 
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
    snpgt = pd.read_csv(SNPGTFILE, sep="\t",low_memory=False)
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
    print expr_annot.shape
    #open output files
    Errorfile = open(TMPDIR+"/Errorfile.out", 'w')
    OUT = open(OUTFILE, "w")
    OUT.write("\t".join(["chrom", "gene", "best.str.start", "best.str.score", "top.variant", "top.variant.score","top.snp.score"])+'\n')
# For each gene, get all cis-variants and the best STR
    for i in range(expr_annot.shape[0]):
        gene=expr_annot.index.values[i]
        ensgene = expr_annot["gene.id"].values[i]  #'ENSG00000215912.7'
        genedir=TMPDIR+"/%s"%gene
        if not os.path.exists(genedir):
            os.mkdir(genedir)
        PROGRESS("Getting data for %s"%gene, printit=DEBUG)
        start = expr_annot["gene.start"].values[i]
        end = expr_annot["gene.stop"].values[i]
    # Pull out cis SNPs
        PROGRESS("Getting cis SNPs for %s"%gene)
        cis_snps = snps[(snps["str.start"] >= (start-DISTFROMGENE)) & (snps["str.start"] <= (end+DISTFROMGENE))]
        print cis_snps.shape , '###'
        cis_snps = cis_snps.loc[cis_snps["str.start"].isin(list(snpgt["start"]))]  ###
        print cis_snps.shape , '###'
        cis_variants = cis_snps.loc[cis_snps['gene']==ensgene]
        cis_variants=cis_variants.sort_values(by="p.wald").head(n=100)
        cis_variants.index = cis_variants["str.start"].apply(lambda x: "SNP_%s"%int(x))
        L=list(cis_variants.index)
        print 'length of cis variants...',len(L) ###############################
    # Pull out most significant STR
        PROGRESS("Getting most significant cis STR for %s"%gene)
        T = strs[strs["gene"]==ensgene].sort_values("p.wald")
        if T.shape[0]==0:
            PROGRESS("There are no STRs found for %s... Gene not in LR table"%gene)
            continue
        else: 
            best_str_start = int(T["str.start"].values[0])
        cis_strs = strs.loc[strs['gene']==ensgene]
        cis_strs.index = list(cis_strs['str.id'])
        try:
            del cis_strs['ID']
        except:
            pass
        cis_variants.loc['STR_'+str(best_str_start)] = list(cis_strs.loc['STR_'+str(best_str_start)])
        print 'best STR ...', '\t', cis_variants.shape
    # Make z file data
        Ztable = MakeZScoreTable(cis_variants[['beta','beta.se']])
        if Ztable is None:
            Errorfile.write(gene+": Z score could not be calculated; beta.se is probably 0 or null\n")
            continue
        else:
            Ztable.to_csv(genedir+'/ZFILE', sep='\t',header=None)
    # Make LD file
        genotypes = snpgt.loc[L]
        genotypes.loc['STR_'+str(best_str_start)] = list(strgt.loc['STR_'+str(best_str_start)])
        del genotypes['chrom']
        del genotypes['start']
        print genotypes.shape
        CorrMatrix = WriteCorrTable(genotypes)
        CorrMatrix.to_csv(genedir+'/LDFILE', sep='\t',header=None, index=None)
        PROGRESS("Matrix of corr was sent to file for %s"%gene)
    #Run caviar
        caviar_cmd = "CAVIAR -l %s -z %s -o %s/caviar -c 1 -f 1 > %s"%(genedir+"/LDFILE", genedir+"/ZFILE", genedir, genedir+"/log")
        os.system(caviar_cmd)

    # Output results
        if not os.path.exists(genedir+'/caviar_post'):
            Errorfile.write(gene+": CAVIAR did not run.\n\tERROR: Segmentation fault (core dumped) in log file\n")
            continue
        else:
            post = pd.read_csv(genedir+'/caviar_post', sep="\t", header=None)
            print 'post','......', post.loc[post[0]=='STR_'+str(best_str_start)][2].tolist()
            caviarstr  =  post.loc[post[0]=='STR_'+str(best_str_start)][2].tolist()[0]
            topvariant =  post.sort_values(post.columns[2], ascending=False).values[0][0]
            topscore  =  post.sort_values(post.columns[2], ascending=False).values[0][2]
            if 'STR' in topvariant:
                snpscore  =  post.sort_values(post.columns[2], ascending=False).values[1][2]
            else:
                snpscore=topscore
            OUT.write("\t".join([CHROM, gene, str(best_str_start), str(caviarstr), topvariant, str(topscore), str(snpscore)])+'\n')
        O=0
    OUT.close()
    Errorfile.close()
