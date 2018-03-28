#!/usr/bin/env python

import argparse
from datetime import datetime
import numpy as np
import scipy.stats
import pandas as pd
import subprocess
import vcf
import sys
import os
import io

"""
    Merge vcf files, format for downstream  and filters homopolymers and seg. dupl STRs. 
    Output is a vcf file to be indexed with tabix
    Usage:
        STR_filter.py --vcf VCFs --homopolymers --seg_dup --hrun --hwe 0.1 --call_rate --heterozygosty 0.3  --out Filtered_STRs.vcf   
"""

SEGDUP="/storage/resources/dbase/human/hg19/hg19_segmentalduplications.bed"
HRUN = "/storage/resources/dbase/human/hg19/hg19.hipstr_reference_hrun.bed"
OUTPUTFILE = ""

HOM = True
SeD = True
HRN = True
HWE = 0.05
ClR = 0.8
HETZYG = 0.3

def PROGRESS(msg):
    sys.stderr.write("%s\n"%msg.strip())

def removehomopolymers(Frame):
    cleanF=Frame.loc[Frame["UNIT"]!=1]
    return(cleanF)

def removeoverlap(Frame, feat ):
    L=list(set(list(Frame['CHROM'])))
    fragments=[]
    t=0
    for C in L:
        X = Frame.loc[Frame['CHROM']==C]
        Y = feat.loc[feat['CHROM']=='chr'+str(C)]
        X['POS'] = X["POS"].astype(int)
        X['END'] = X["END"].astype(int)
        for i in range(len(list(Y.index))):
            start = list(Y['START'])[i]
            end = list(Y['END'])[i]
            X2 = X.loc[(X["END"]<=start) | (X["POS"]>=end)]
            X = 0; X = X2
        fragments.append(X2.sort_values('POS'))
        print(C,'\t',X.shape)
    result = pd.concat(fragments)
    return(result)
   
def removelowcallrate(Frame): 
    Frame['Count0'] = Frame.isnull().sum(axis=1)
    Frame['Count1'] = Frame.isin({'./.:.'}).sum(1)
    Frame['New'] = 650 - Frame['Count0']                     #650 samples
    result = Frame.loc[Frame['Count1']<Frame['New']*0.2]     #Call rate 80%
    del result['Count0']
    del result['Count1']
    del Table['New']
    retun (result)
    
def GetLocusStats(record, samples=[]):
    hwe_p = 0
    het = 0
    # Get genotypes, allele frequencies
    allele_counts = {}
    obs_het = 0
    obs_hom = 0
    total = 0
    for sample in record:
        if len(samples)>0 and sample.sample not in samples: continue
        if sample["GB"] == "." or sample["GB"] == None: continue
        gt = map(int, sample["GB"].split("|"))
        if gt[0] == gt[1]: obs_hom += 1
        else:
            obs_het += 1
        total += 1
        for al in gt:
            allele_counts[al] = allele_counts.get(al, 0) + 1
    # Get Allele frequencies
    allele_freqs = {}
    for key in allele_counts.keys():
        allele_freqs[key] = allele_counts[key]*1.0/sum(allele_counts.values())
    # Get expected num homs/hets
    exp_hom_frac = 0
    for al in allele_freqs.keys():
        exp_hom_frac += allele_freqs[al]**2
    # Binomial test for HWE
    hwe_p = scipy.stats.binom_test(obs_het, n=obs_het+obs_hom, p=1-exp_hom_frac)
    # Compute heterozygosity
    het = 1-sum([allele_freqs[al]**2 for al in allele_freqs.keys()])
    # Get mean allele length
    mean_allele = sum([al*allele_freqs[al] for al in allele_freqs])
    return (hwe_p, het, mean_allele,obs_het+obs_hom)

def getper(a, infofield):
    a = a.split(";")
    b =[b.split("=")[1] for b in a if infofield in b]
    return(b[0])

def addheader():
    header = "\n".join(['##INFO=<ID=HWE,Number=1,Type=Float,Description="HWE pvalue genotype frequencies not as expected">',
            '##INFO=<ID=HET,Number=1,Type=Float,Description="Heterozygosity">',
            '##INFO=<ID=CCOUNT,Number=1,Type=Float,Description="Number of samples with genotype information">',
            '##FILTER=<ID=HET,Description="Heterozygosity less than '+str(HETZYG)+'">',
            '##FILTER=<ID=HRUN,Description="Hrun greater than -1">',
            '##FILTER=<ID=HWE,Description="HWE less than '+str(HWE)+'">',
            '##FILTER=<ID=CALLRATE,Description="Callrate less than '+str(ClR)+'">',
            '##FILTER=<ID=HOM_POLY,Description="Homopolymer locus">',
            '##FILTER=<ID=SEGDUP,Description="Locus in a segmental duplication">\n#'])
    return(header)

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF file", type=str, required=True)
    parser.add_argument("--unrelated-samples", help="Restrict to these samples to calculate popgen stats", type=str, required=False)
    parser.add_argument("--out", help="Write filtered data to this file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict to chromosome", type=str, required=False)
    parser.add_argument("--debug", help="Print debug info", action="store_true")
    parser.add_argument("--homopolymers", help="Filters and Remove homopolymers STRs", action="store_true")
    parser.add_argument("--seg_dup", help="STRs overlapping segmental duplication will be removed", action="store_true")
    parser.add_argument("--hrun", help="Hexa and penta STRs with long homopolymer runs will be removed", action="store_true")
    parser.add_argument("--hwe", help="Remove locus not at hwe is pval >= 0.05 or otherwise specified ", type=float, required=False)
    parser.add_argument("--call_rate", help="Remove locus with call rate <= 80%% or otherwise specified ", type=float, required=False)
    parser.add_argument("--heterozygosity", help="Remove locus with low heterozygosity <=0.3 or otherwise specified ", type=float, required=False)

#Input variables
    args = parser.parse_args()
    VCF = args.vcf
    if args.chrom:
        if 'chr' in args.chrom:
            CH=args.chrom[3:]
        else:
            CH=args.chrom
    HOM = args.homopolymers
    SeD = args.seg_dup
    HRN = args.hrun
    HWE = args.hwe
    ClR = args.call_rate
    if args.heterozygosity:
        try:
            HETZYG = float(args.heterozygosity)
        except:
            HETZYG = 0.3
    if args.unrelated_samples:
        samples = [item.strip() for item in open(args.unrelated_samples, "r").readlines()]
    else:
        sample = []
    OUTPUTFILE = args.out
    
#Setup
    PROGRESS('Starting ... ')
    with open(VCF, 'r') as f:
        lines = [l for l in f if not l.startswith('##')] 
        Head = [l for l in f if l.startswith('##')]
    with open(OUTPUTFILE, 'w') as f:
        f.write('\n'.join(Head))
               
    Table = pd.read_table( io.BytesIO(str.join(os.linesep, lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}).rename(columns={'#CHROM': 'CHROM'})
        
    vcf_reader = vcf.Reader(filename=VCF)
    N=len(vcf_reader.samples)*ClR
    
    if args.chrom:
        Table = Table.loc[Table['CHROM']==CH]
        vcf_reader = vcf.Reader(filename=VCF).fetch(CH)
    PROGRESS('File opened... '+str(Table.shape))
    
    Table['REF'] = [x.split(',')[0] for x in list(Table['REF'])]
    Table['UNIT'] = [int(getper(x, 'PERIOD')) for x in list(Table['INFO']) ]
    Table['END'] = [int(getper(x,'END')) for x in list(Table['INFO'])]    
    PROGRESS("unit & end"+str(Table.shape))
    
    Table["stats"]=[GetLocusStats(record) for record in vcf_reader]
    Table["hwepval"] = [x[0] for x in list(Table['stats'])]
    Table["Hetzyg"] = [x[1] for x in list(Table['stats'])]
    Table["mean_al"] = [x[2] for x in list(Table['stats'])]
    Table["CallCounts"]=[x[3] for x in list(Table['stats'])]               
    PROGRESS('Additional details added... '+str(Table.shape[0]))
    
# low call rate STRs and homopolymers
# low heterozygosity locus            
# locus not in HWE


    if HOM:
        Table["FILTER"] = np.where(Table["UNIT"] !=1, Table["FILTER"], Table["FILTER"]+" HOM_POLY")
    ### Create HOM only file with ==1 instead of !=1

    if ClR:     
        Table["FILTER"] = np.where(Table["CallCounts"] >= N, "", "CALLRATE")
    
    if HWE:
        Table["FILTER"] = np.where(Table["hwepval"] >=HWE, Table["FILTER"], Table["FILTER"]+" HWE")
        
    if HETZYG:
        Table["FILTER"] = np.where(Table["Hetzyg"] >=HETZYG, Table["FILTER"], Table["FILTER"]+" HET")
        
    PROGRESS('Call rate, homopolymers, hwe, heterozyg. Filters added... '+str(Table.shape[0]))    
    
#Seg dup
    if SeD:
        Seg_dup = pd.read_csv(SEGDUP, sep='\t', header=None)
        Seg_dup.columns = ['CHROM', 'START','END','OTHERS','INFO','STRAND']
        Table_c = removeoverlap(Table, Seg_dup)
        Table["FILTER"] = np.where(Table["ID"].isin(list(Table_c['ID'])), Table["FILTER"], Table["FILTER"]+" SEGDUP")    
    PROGRESS('Segmental Duplication filter added... '+str(Table.shape))
# Hrun
    if HRN:
        X = pd.read_csv(HRUN,sep='\t', header=None)
        X.columns = ['CHROM', 'START','END','Unitsize','maxrun']
        X1 = X.loc[X['Unitsize'].isin([5,6])]
        hrun = X1.loc[X1['maxrun']>X1['Unitsize']]
        Table_c = removeoverlap(Table, hrun)
        Table["FILTER"] = np.where(Table["ID"].isin(list(Table_c['ID'])), Table["FILTER"], Table["FILTER"]+" HRUN ")
    PROGRESS('Penta and hexa with homopolymer runs ... '+str(Table.shape))

    
#Reformat INFO, FILTER    
    Table["info2"] = ';HWE='+Table["hwepval"].astype(str) +';HET=' +Table["Hetzyg"].astype(str)+ ';CCOUNT='+ Table["CallCounts"].astype(str)  
    Table['INFO'] = Table['INFO']+Table['info2']
    Table['FILTER'] = np.where(Table["FILTER"]!='', Table["FILTER"], "PASS")
    C = [';'.join(s.split()) for s in list(Table['FILTER'])]
    Table['FILTER'] = C
    PROGRESS('Fields formatted ... ')
    
#Clean up    
    del Table['stats']
    del Table['hwepval']
    del Table['mean_al']
    del Table['Hetzyg']
    del Table['CallCounts']
    del Table['UNIT']
    del Table['END']
    del Table['info2']
    C = ['-'.join(x.split('-')[:2]) for x in list(Table.columns)]
    Table.columns = C
    C = Table.loc[Table['FILTER']=='PASS']
    PROGRESS('Cleaned up. Saving to file ... \n\t'+str(C.shape[0])+' passed all filters...')

#Package and save vcf
    PROGRESS("Saving to file")
    command = "grep '^##' "+VCF
    vcfheader = subprocess.check_output(command, shell=True)
    f=open('tmp','w')
    f.write(vcfheader.decode('utf-8'))
    f.write(addheader())
    f.close()
    Table = Table.sort_values(['CHROM','POS'])
    Table.to_csv('table.tab',sep='\t',index=None)
    command = "cat tmp table.tab >"+OUTPUTFILE 
    MG = subprocess.check_output(command, shell=True)
    command = "rm table.tab"
    MG = subprocess.check_output(command, shell=True)
    command = "rm tmp"
    MG = subprocess.check_output(command, shell=True)
#compress and index vcf 
    PROGRESS("Indexing")
    command = "bgzip -c "+ OUTPUTFILE +" > "+OUTPUTFILE+'.gz'
    MG = subprocess.check_output(command, shell=True)
    command = "tabix -p vcf "+ OUTPUTFILE +'.gz '
    MG = subprocess.check_output(command, shell=True)


    
if __name__ == "__main__":
    main()

