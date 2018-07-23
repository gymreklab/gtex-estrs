#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import sys

"""
Merge causality scores from all chromosome
Usage: ./path/to/merge_caviar_scores.py <DATADIR> <tissuelist>
Output will be put in DATADIR as "CAVIAR_wg.csv" and gene_str_score.tab file with all scores calculated for each gene
Prefix of all intermediate files is "caviar.out  followed by chrom number
"""

try:
    DATADIR = sys.argv[1]      #Tissues Directory. - Assumption is caviar intermidiate files are in HH DIR within this DIR
    TISSUES = sys.argv[2].split(",") #Tissue(s)
except:
    sys.stderr.write(__doc__)
    sys.exit(1)


#Merge all chromosome and move intermediate files into a new dir       ##Not done 2,3
chrom=[i for i in range(1,23,1)]

for T in TISSUES:
    path='%s/%s/HH/caviar.out'%(DATADIR,T)
    frames = []
    scores = []
    
    command='mv %s/%s/HH/CAVIAR_analysis.table %s/%s/HH/old_CAVIAR.table'%(DATADIR, T,DATADIR,T) #can remove this
    os.system(command)
    
    for C in chrom:
        print(C)
        frames.append(pd.read_csv('%s%s'%(path,str(C)), sep='\t'))
        scores.append(pd.read_csv('%s/%s/HH/strs_scorechr%s'%(DATADIR, T,str(C)), sep='\t'))
        
    data = pd.concat(frames)
    data.to_csv('%s/%s/HH/CAVIAR_wg.csv'%(DATADIR,T), sep='\t', index=None)
    #
    SCOR = pd.concat(scores)
    SCOR.to_csv('%s/%s/HH/gene_str_score.tab'%(DATADIR,T), sep='\t', index=None)
    
    command='mkdir %s/%s/HH/intermediate_files_caviar'%(DATADIR, T)
    os.system(command)
    
    command='mv %s/%s/HH/caviar.out* %s/%s/HH/intermediate_files_caviar/'%(DATADIR, T,DATADIR,T)
    os.system(command)
    command='mv %s/%s/HH/strs_scorechr* %s/%s/HH/intermediate_files_caviar/'%(DATADIR, T,DATADIR,T)
    os.system(command)
    print(T, SCOR.shape)
print("END")

#'/home/szfeupe/projects/GTEX_eSTRs/gtex-estrs/Scripts/STRs-vs-SNPs_effects/merge_caviar_scores.py'