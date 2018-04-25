#!/usr/bin/env python2.7

# Libraries
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vcf

"""
Merging the VCF samples.
This script Assume the same INFO column for all samples which is false
"""

VCF='/storage/resources/datasets/gtex/vcfs_650/'
SAMPLES='/storage/resources/datasets/gtex/GTExRunTable.txt'   ##column 33 has sample IDs. fROM this we have hibrid PCR and PCR free

#Assessment of sample size
runtable = pd.read_csv(SAMPLES, sep='\t', header=0)
data = runtable[['Assay_Type','SRA_Sample', 'Run','Sample_Name','body_site', 'data_type','histological_type','analyte_type']]
data['SampleID']=data.Sample_Name.apply(lambda x: '-'.join([x.split('-')[0],x.split('-')[1],x.split('-')[2]]))
data_rna = data.loc[data['Assay_Type']=='RNA-Seq']
data_wgs = data.loc[data['Assay_Type']=='WGS']

samplesrun = os.listdir(VCF)
vcfs = [x for x in samplesrun if 'hipstr.vcf.gz' in x and '.tbi' not in x]
run1 = [x.split('_')[0] for x in vcfs]
wgs1 = data_wgs.loc[data_wgs['Run'].isin(run1)]
print 'There are ',len(vcfs),' wgs in DIR to be processed'

##  Merging of the genotypes into one file
i=0
# First file sets the tone
    #  Header
command=["zgrep '^##' "+VCF+vcfs[0] + " > subset_strs.vcf"]
output=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

    #  INFO column
df = pd.read_csv(VCF+vcfs[0], compression='gzip', header=0, skiprows=int(123), sep='\t', quotechar='"', low_memory=False)
df['#CHROM']=df['#CHROM'].astype(str)
STRs = df[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]
print 'First sample IN'

    # Add up samples
for files in vcfs[1:]:
    command=["zgrep '^##' "+VCF+files + " |wc -l"]
    output=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (N, err) = output.communicate()
    df = pd.read_csv(VCF+files, compression='gzip', header=0, skiprows=int(N), sep='\t', quotechar='"', low_memory=False)
    df['#CHROM']=df['#CHROM'].astype(str)
    if len(set(list(df['#CHROM'])))<24:
        print files
        continue
    del df['INFO']
    del df['QUAL']
    del df['FILTER']
    STRs = pd.merge(STRs,df,how='left', on=['#CHROM','POS','ID','REF','ALT','FORMAT'])
    i+=1

STRs.to_csv("/storage/szfeupe/Runs/650GTEx_estr/Genotypes/STRs", sep='\t',index=None)

#command=["cat ~/projects/GTEX_eSTRs/650GTEX/Preprocessing_Scripts/subset_strs.vcf STRs > /storage/szfeupe/Runs/650GTEx_estr/Genotypes/AllSamplesSTRs.vcf"]