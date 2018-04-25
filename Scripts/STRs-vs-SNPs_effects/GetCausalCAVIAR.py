import argparse
import math
import numpy as np
import os
import pandas as pd
import random
import shutil
import sys
import gzip

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
    for V1 in variants:
        COV=[]
        for V2 in variants:
            X=G[V1].replace('None', np.nan).astype(float)
            Y=G[V2].replace('None', np.nan).astype(float)
            if X.corr(Y) is np.nan:    #### 
                COV.append(0.0) #For missing LD we assume non linear corr (undetermined LD)
            else:
                COV.append(X.corr(Y))
        CMat.append(COV)
    return pd.DataFrame(CMat,columns=variants, index=variants) 

def lookfor (x,p):
    for row in range(1,len(p.index)):
        if x in p.values[row][0]:
            top = p.values[row][0]
            score = p.values[row][2]
            return(row,top, score)


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
    DEBUG =False
    PROGRESS("\nLoad expression", printit=DEBUG)
    expr = pd.read_csv(EXPRFILE)
    samples_to_keep = list(expr.index)
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
    del snps['Unnamed: 0']
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
    str_samples = samples_to_keep
    expr = expr.reindex(str_samples)
    snpgt = snpgt[["chrom","start"] + str_samples]
    snpgt.index = list(snpgt["start"].apply(lambda x: "SNP_%s"%int(x)))
    strgt = strgt[["chrom","start"] + str_samples]
    strgt.index = list(strgt["start"].apply(lambda x: "STR_%s"%int(x)))
    # Load eSTR results
    PROGRESS("Restrict to eSTR genes only", printit=DEBUG)
    if ESTRGENESFILE is not None:
        estr_genes = pd.read_csv(ESTRGENESFILE, sep="\t")
        Genes = estr_genes.loc[estr_genes['qvalue']<=0.1]['gene']  # estrs at 10%FDR
        expr_annot = expr_annot.loc[expr_annot['gene.id'].isin(list(Genes))]
    #open output files
    Errorfile = open(TMPDIR+"/Errorfile.out", 'w')
    OUT = open(OUTFILE, "w")
    PROGRESS("Start output file "+OUTFILE, printit=DEBUG)
    OUT.write("\t".join( ['CHROM','gene','num.strs.in.top5','top_snp','top_snp_score','top_str','top.str.score','str.rank'])+'\n')
    # For each gene, get all cis-variants and the best STR
    O=0
    for i in range(expr_annot.shape[0]):
        gene=expr_annot.index.values[i]
        ensgene = expr_annot["gene.id"].values[i]  #'ENSG00000215912.7'
        genedir=TMPDIR+"/%s"%gene
        if not os.path.exists(genedir):
            os.mkdir(genedir)
        clear_cmd = "rm "+genedir+'/*'
        os.system(clear_cmd)
        PROGRESS("Getting data for %s"%gene, printit=DEBUG)
        start = expr_annot["gene.start"].values[i]
        end = expr_annot["gene.stop"].values[i]
    # Pull out cis SNPs
        PROGRESS("Getting cis SNPs for %s"%gene)
        cis_snps = snps[(snps["str.start"] >= (start-DISTFROMGENE)) & (snps["str.start"] <= (end+DISTFROMGENE))]
        cis_snps = cis_snps.loc[cis_snps['gene']==ensgene]
        cis_snps.index = cis_snps["str.start"].apply(lambda x: "SNP_%s"%int(x))
        print (cis_snps.shape , '##SNPs#')
        cis_variants = cis_snps.loc[cis_snps["str.start"].isin(list(snpgt["start"]))]  ###        
        #cis_variants = cis_snps.loc[cis_snps['gene']==ensgene]
        cis_snps=cis_variants.sort_values(by="p.wald").head(n=100)
        #cis_snps.index = cis_snps["str.start"].apply(lambda x: "SNP_%s"%int(x))
        L=list(cis_snps.index)
    # Pull out cis STR

    # Pull out cis STR
        PROGRESS("Getting most significant cis STR for %s"%gene)
        cis_strs = strs[strs["gene"]==ensgene].sort_values("p.wald")
        if cis_strs.shape[0]==0 :
            PROGRESS("There are no STRs found for %s... Gene not in LR table"%gene)
            continue
        elif cis_snps.shape[0]<=1:
            PROGRESS("There are no or not enough SNPs found for %s... Gene not in LR table"%gene)
            continue
        else: 
            cis_strs.index = cis_strs["str.start"].apply(lambda x: "STR_%s"%int(x))
            best_str_start = int(cis_strs["str.start"].values[0])
            L0 = list(cis_strs.index)
        #
        cis_variants = pd.concat([cis_snps, cis_strs])
        print(len(L), len(L0), cis_variants.shape)
    # Make z file data
        Ztable = MakeZScoreTable(cis_variants[['beta','beta.se']])
        if Ztable is None:
            Errorfile.write(gene+": Z score could not be calculated; beta.se is probably 0 or null\n")
            continue
        else:
            Ztable.to_csv(genedir+'/ZFILE', sep='\t',header=None)
    # Make LD file
        genotypes = snpgt.loc[L]
        genotypes = pd.concat([genotypes, strgt.loc[L0] ])
        del genotypes['chrom']
        del genotypes['start']
        CorrMatrix = WriteCorrTable(genotypes)
        CorrMatrix.to_csv(genedir+'/LDFILE', sep='\t',header=None, index=None)
        PROGRESS("Matrix of corr was sent to file for %s"%gene)        
    #Run caviar
        caviar_cmd = "CAVIAR -l %s -z %s -o %s/caviar -c 1 -f 1 > %s"%(genedir+"/LDFILE", genedir+"/ZFILE", genedir, genedir+"/log")
        os.system(caviar_cmd)
    #Output results
        if not os.path.exists(genedir+'/caviar_post'):
            Errorfile.write(gene+": CAVIAR did not run.\n\tERROR: Segmentation fault (core dumped) in log file\n")
            continue
        else:
            post = pd.read_csv(genedir+'/caviar_post', sep="\t", header=None)
            post = post.sort_values(post.columns[2], ascending=False)
            post = post.reset_index(drop=True)
            p = post.head(5)
            #print (p[0])
            num_str = len([x for x in list(p[0].values) if "STR_" in x])
            if 'STR_' in p[0][0]:
                topstr = p[0][0]
                topstrscore = p.values[0][2]
                I,topsnp, topsnpscore =lookfor('SNP_', post)
            else:
                topsnp = p[0][0]
                topsnpscore = p.values[0][2]
                I, topstr , topstrscore =lookfor('STR_',post)
            OUT.write("\t".join([CHROM, gene, str(num_str),topsnp,str(topsnpscore),str(topstr), str(topstrscore),str(I+1)])+'\n')
            strsscores = post.loc[post[0].isin(cis_strs.index)][[0,2]]
            strsscores['chrom']=[CHROM]*strsscores.shape[0]
            strsscores['gene']= [gene]*strsscores.shape[0]
            strsscores.columns = ['str','score', 'chrom', 'gene']
            with open('strs_score'+CHROM, 'a') as k:
                (strsscores[['chrom', 'gene','str','score']]).to_csv(k, header=False, index=False, sep='\t')
    #   #
            #
    OUT.close()
    k.close()
    Errorfile.close()
    print O
