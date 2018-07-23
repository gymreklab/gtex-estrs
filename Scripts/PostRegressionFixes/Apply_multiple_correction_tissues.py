#! /bin/python2.7

import pandas as pd
import numpy as np
import subprocess
import collections
import os

## FDR correction After gene level adjustment
SHORTEN = {
    "Artery-Aorta":"Artery A."     ,
    "Artery-Tibial": "Artery T.",
    "Adipose-Subcutaneous": "Adipose S.",    
    "Adipose-Visceral":"Adipose V.",
    "Brain-Caudate":"Caudate"   , 
    "Brain-Cerebellum":"Cerebellum",
    "Cells-Transformedfibroblasts": "Fibroblast",
    "Esophagus-Mucosa": "E. Mucosa",
    "Esophagus-Muscularis":"E Muscularis",
    "Heart-LeftVentricle":"Ventricule",
    "Lung": "Lung",
    "Muscle-Skeletal": "Muscle",
    "Nerve-Tibial":"Nerve",
    "Skin-NotSunExposed": "Skin Unexposed",
    "Skin-SunExposed":"Skin Leg",
    "Thyroid":"Thyroid",
    "WholeBlood": "Blood"
}

tissu = sorted(SHORTEN.keys())
#path = "/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/"
path = "/storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/"

def fdrcorrection(tissue):
    print tissue, ' variants ...'
#Get most signif. variant by gene from linear reg STRs
    LR1=pd.read_csv(path+ tissue+"/Lin_Reg_Out" , '\t')
    
#Locus level
    LR1['p.wald'].to_csv('pvalues.txt', sep='\n', index=False)
    Tell = subprocess.call("/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/PostRegressionFixes/fdr-correct.r")
    Qval=pd.read_csv('/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/PostRegressionFixes/qvalues.txt', sep=' ')
    
    LR1['llqvalue']=list(Qval['qvalue'])
    LR1['llsignif']=list(Qval['significant'])

#Gene level
    LR0 = LR1.sort_values("p.wald").groupby("gene", as_index=False).first()     
    print(LR1.shape, '  to  ', LR0.shape)

    #Add counts tests by gene
    counts=pd.DataFrame({'cts' : LR1.groupby(["gene"]).size()})    ## This is the count by genes
    genes = list(LR0['gene'])
    LR0['NTEST']= list(counts.loc[genes]['cts'])
    
    #Gene level adjustment
    #(1) min_pval* #test
    LR0['AD.pval']=LR0['p.wald']*LR0['NTEST']
    #(2) if AD_pval>1 => AD_pval=1
    LR0['AD.pval'][LR0['AD.pval']>1] = 1
    
    #Save pval in file and FDR correct
    LR0['AD.pval'].to_csv('pvalues.txt', sep='\n', index=False)
    Tell = subprocess.call("/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/PostRegressionFixes/fdr-correct.r")
        
    #FDR corrected... add to dataframe
    Qval=pd.read_csv('/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/PostRegressionFixes/qvalues.txt', sep=' ')
    LR0['qvalue']=list(Qval['qvalue'])
    LR0['significant']=list(Qval['significant'])

#Merging
    merging=['gene','chrom','str.id','str.start','beta','beta.se','p.wald','llqvalue','llsignif']
    LRP = pd.merge(LR1,LR0, on=merging, how='left')
    
#Header arrangement
    Head=['gene','chrom','str.id','str.start','p.wald','llqvalue','llsignif','NTEST','qvalue','significant','beta','beta.se']
    Out=LRP[Head]
    Out.to_csv(path+tissue+'/PQValues', sep='\t', index=False)

    S=LR0['AD.pval']
    print len(S),' total tests... ', len(S[S>=1]) , ' pvalues were reduced to 1'
    print len(LRP[LRP['qvalue'] <=0.1]),'\t gene level qval<=0.1')
    print len(LRP[LRP['llqvalue'] <=0.1]),'\t locus level qval<=0.1'
    print len(LRP[LRP['llqvalue'] <=0.01]),'\t qval<0.01\n'
    return()
#
#
for T in tissu[1:]:                     
    fdrcorrection(T+'/SNP_Analysis')
    #fdrcorrection(T+'/SNP_Analysis')