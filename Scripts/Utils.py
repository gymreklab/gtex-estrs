#!/usr/bin/env python2.7
import matplotlib.pyplot as plt
import pandas as pd


"""
                UTIL1
Merge linear regression output from all chromosomes
USAGE:
    mergeallchrs('Cells-Transformedfibroblasts/SNP_Analysis')
"""

def mergeallchrs(tissue):
    Frames=[pd.read_csv('/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+ tissue+'/chr1/Lin_Reg_Out',sep='\t')]
    for x in range(2,23):
        LN='/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+ tissue+'/chr'+str(x)+'/Lin_Reg_Out'
        frame1=pd.read_csv(LN, sep='\t')
        Frames.append(frame1)
#print('Chr',x,'    ',frame1.shape)
    Results=pd.concat(Frames)
    print('\n All Chrms ','   ', Results.shape)
    Results.to_csv('/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+ tissue+'/Lin_Reg_Out', sep='\t', header=True)
    print('Continuing...')

    Frames=[pd.read_csv('/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+tissue+'/chr1/Lin_Reg_Out_perm',sep='\t')]
    for x in range(2,23):
        LN='/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+ tissue+'/chr'+str(x)+'/Lin_Reg_Out_perm'
        frame1=pd.read_csv(LN, sep='\t')
        Frames.append(frame1)
#print('Chr',x,'    ',frame1.shape)
    Results=pd.concat(Frames)
    print('\n All Chrms ','   ', Results.shape)
    Results.to_csv('/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+ tissue+'/Lin_Reg_Out_perm', sep='\t', header=True)

    print('THE END')
    

