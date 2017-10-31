import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import subprocess

def round1adjust(annot, linreg):
    #Set up empty dataframe
    index=['0']
    Out_Lin_reg = pd.DataFrame(index=index, columns=linreg.columns)
    Out_Lin_reg = Out_Lin_reg.fillna(0)
    Out_Lin_reg['PVAL']=0
    Test=[]        #Number of STRs tested for a gene
    Notest=[]      #Gene twith no tested genes
    for index, gene in annot.iterrows():
        geneid = gene['gene.id']
        start = gene['gene.start']
        stop = gene['gene.stop']
#        strs = linreg.loc[linreg['str.start'].isin(range(start,stop))]
        strs = linreg.loc[linreg['gene'].isin([geneid])]
        N = len(strs)
#identifying the STR test with smallest pvalue
        if N>0 :
            ind_low = strs.loc[strs['p.wald'].idxmin()]
            lowpval = strs.loc[[ind_low.name]]
#adjusting by the number of tests
            Adj_pval = lowpval['p.wald']*N
            lowpval['PVAL'] = Adj_pval 
#Append product in the output for final FDR adjustment
            Out_Lin_reg.loc[len(Out_Lin_reg)] = lowpval.values.tolist()[0]
            Test.append(N)
    
        else:
            Notest.append(gene['gene.id']) 
    Out_Lin_reg = Out_Lin_reg.drop('0')
    print(len(Notest), ' genes were not tested for eSTRs')
    print(len(Test), ' genes were tested for eSTRs\n') 
    return(Out_Lin_reg, Notest, Test)


linregout=str(sys.argv)[1]
Input='/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'

#open gene annotation  gene.chr	gene.start	gene.stop	gene.id
Annot = pd.read_csv('~/projects/GTEX_eSTRs/data/Lin_Reg/Gene_Exp_Annotation.txt', sep=',')
#open Linear regression output    gene	chrom	str.start	p.wald
Linreg = pd.read_csv(linregout, sep='\t')
#Chromosome list set
chrom = ['chr'+str(i) for i in range(1,22)]
chrom.append('chrX') ; chrom.append('chrY')
#Setting up output
index=['0']
Adjusted = pd.DataFrame(index=index, columns=Linreg.columns)
Adjusted = Adjusted.fillna(0)
print (Adjusted.shape)
NT = []


#    Single out the STR with lowest pval and FDR adjust by chromosome
for ch in chrom:
    A = Annot.loc[Annot['gene.chr'].isin([ch])]
    LR = Linreg.loc[Linreg['chrom'].isin([ch])]
    print(ch,' ',len(A), len(LR))
    Adjval, GTEST, NTest = round1adjust(A, LR)
#    Adjval.drop('0')
    Adjusted = pd.concat([Adjusted, Adjval])#Adjusted.append(Adjval, ignore_index=True)
    NT =NT + NTest
    print('Done with ', ch, '\t', Adjval.shape,' TO... ',Adjusted.shape)
print('End')

#    Save Linear reg. + FDR adjusted in output
Adjusted.to_csv(linregout, sep='\t')

print(len(NT))

Out_Lin_reg=Adjusted
Out_Lin_reg=Out_Lin_reg.drop('0')
#SOmetime, the p-values multiplied by the number of STRs is >1 
Out_Lin_reg['NTest']=NT
Out_Lin_reg['PVAL'][Out_Lin_reg['PVAL']>1] = 1
PVAL=Out_Lin_reg['PVAL']
PVAL.to_csv('pvalues.txt', sep='\n', index=False)

print(len(PVAL))
print(len(PVAL[PVAL>=1]))

#Now, we use the QVALUE package to adjust the pvalues and obtain the qvalues
Tell = subprocess.call("/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/fdr-correct.r")
Qval=pd.read_csv('/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/qvalues.txt', sep=' ')
print (Qval.shape, '\t', Out_Lin_reg.shape)
Out_Lin_reg['pvalue1']=list(Qval['pvalue'])
#Out_Lin_reg['#oftest']=Test
Out_Lin_reg['qvalue']=list(Qval['qvalue'])
Out_Lin_reg['significant']=list(Qval['significant'])

Out_Lin_reg.to_csv(Input+'PQValues.txt', sep='\t', index=False)
#print(Out_Lin_reg[Out_Lin_reg['pvalue'] >=1].shape)
print("FDR correction summary: \neSTRs counts\t Treshold")
print(len(Qval[Qval['qvalue'] <=0.1]),'\t qval<=0.1')
print(len(Qval[Qval['qvalue'] <=0.05]),'\t qval<=0.05')
print(len(Qval[Qval['qvalue'] <=0.01]),'\t qval<0.01')

print('FDR correction done. Continue with Anova test for Heritability')