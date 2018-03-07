#!/usr/bin/env python2.7
import argparse
import pandas as pd
PATH = "/storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/"

"""
                UTIL
Merge linear regression output from all chromosomes
USAGE:
    mergeallchrs('Cells-Transformedfibroblasts/SNP_Analysis')
"""

def mergeallchrs(tissue, PATH):
    print tissue, '... START'
    Frames=[pd.read_csv(PATH+ tissue+'/chr1/Lin_Reg_Out',sep='\t')]
    print PATH+ tissue+'/chr1/Lin_Reg_Out'
    for x in range(2,23):
        LN=PATH+ tissue+'/chr'+str(x)+'/Lin_Reg_Out'
        print LN
        frame1=pd.read_csv(LN, sep='\t')
        Frames.append(frame1)
#print('Chr',x,'    ',frame1.shape)
    Results=pd.concat(Frames)
    print '\n All Chrms ','   ', Results.shape, PATH+ tissue+'/Lin_Reg_Out'
    Results.to_csv(PATH+ tissue+'/Lin_Reg_Out', sep='\t', header=True)
    print 'Continuing...'

    Frames=[pd.read_csv(PATH+tissue+'/chr1/Lin_Reg_Out_perm',sep='\t')]
    for x in range(2,23):
        LN=PATH+ tissue+'/chr'+str(x)+'/Lin_Reg_Out_perm'
        frame1=pd.read_csv(LN, sep='\t')
        Frames.append(frame1)
#print('Chr',x,'    ',frame1.shape)
    Results=pd.concat(Frames)
    print '\n All Chrms ','   ', Results.shape
    Results.to_csv(PATH+ tissue+'/Lin_Reg_Out_perm', sep='\t', header=True)
    print 'THE END'
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge regression results from all chromosomes")
    parser.add_argument("--dir", help="Analysis directory", type=str, required=True)
    parser.add_argument("--tissue", help="tissue name as it appears in the analysis directory", type=str, required=True)
    args = parser.parse_args()
    PATH = args.dir +'/'   
    tissue = args.tissue
    
    mergeallchrs(tissue, PATH)

