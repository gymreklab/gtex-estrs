import argparse
import math
import numpy as np
import os
import io
import sys
import subprocess
import pandas as pd
import vcf
from datetime import datetime
"""
    Merge vcf files, format for downstream  and filters homopolymers and seg. dupl STRs. 
    Output is a vcf file to be indexed with tabix
    Usage:
        MergeSTRvcf.py --vcfdir DIR_withVCFs --remove_homopolymers --remove_seg_dup --out Merged_STRs.vcf
        OR
        MergeSTRvcf.py --vcfdir DIR_withVCFs --vcflist file/w.listofvcfs --out Merged_STRs.vcf --remove_homopolymers --remove_seg_dup
        The first command will take all files in  DIR_withVCFs
"""
    
VCFDIR = None
OUTFILE = ""
DEBUG = False
tmp = 'temp_files/'
OUT = 'Genotypes/'


def PROGRESS(msg):
    sys.stderr.write("%s\n"%msg.strip())

def combine(X):
    X['LIST']=X.astype(str).apply(lambda x: ','.join(x), axis=1)
    X['L']=X['LIST'].astype(str).apply(lambda x: ','.join([x for x in set(x.strip().split(',')) if str(x)!= 'nan']))
    return(list(X['L']))

def Mergesamples(X):
    PROGRESS('... Merge   %s'%str(datetime.now().strftime('%H:%M:%S')))
    S = list(X.keys())
    MG = X[S[0]][["#CHROM","POS","ID","REF","ALT","QUAL","FILTER", "FORMAT","U.size","END"]]
    for s in S:
        MG = pd.merge(MG, X[s], how='outer', on=["#CHROM","POS","ID","ALT","QUAL","FILTER", "FORMAT","U.size","END"])
    ref = [x for x in list(MG.columns) if x[:3]=='REF']; ref=list(set(ref))
    n=len(ref)
    MG['combined_ref']=combine(MG[ref]) 
    MG['combined_info']=combine(MG[list(set([x for x in list(MG.columns) if x[:4]=='INFO']))])
    for x in ref:
        del MG[x]
    return(MG)

def Getgb(record):
    spl = record.samples
    if spl[0]['GB'] is None:
        gb='.'
    else:
        gb=spl[0]['GB']
    geno = ':'.join(['./.', gb])
    return(geno)

def removehomopolymers(Frame):
    cleanF=Frame.loc[Frame["U.size"]!=1]
    return(cleanF)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get specific arguments for parsing and reformat")
    parser.add_argument("--vcfdir", help="Path Directory containing the vcf files", type=str, required=True)
    parser.add_argument("--vcflist", help="file containing the list of vcfs file names to be merged.    NO PATH. This is helpful if splitting the sample and run in parallel. Did not want to go into multicore", type=str, required=False)
    parser.add_argument("--out", help="Write merged data to this file", type=str, required=True)
    parser.add_argument("--debug", help="Print debug info", action="store_true")
    parser.add_argument("--remove_homopolymers", help="Remove homopolymers STRs", action="store_true")
#    parser.add_argument("--remove_seg_dup", help="STRs overlapping segmental duplication will be removed", action="store_true")

    args = parser.parse_args()

    PROGRESS(".\nAssessing input variables")
    OUTFILE = args.out
    VCFDIR = args.vcfdir
    if args.debug: DEBUG = True
    if args.remove_homopolymers:
        HOMO_POLY=True
    else:
        HOMO_POLY=False
    if args.vcflist is not None:
        VCFLIST = args.vcflist
        vcfs = open(VCFLIST,'r').readlines()
        vcfs = [VCFDIR+s.strip() for s in vcfs]
    else:  
        VCFLIST = os.listdir(VCFDIR)
        vcfs = [VCFDIR+s.strip() for s in VCFLIST]

    attrib = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER", "INFO","FORMAT"]
    CHR =[str(i) for x in range(1,23,1)] +['X','Y']
     
    # Merging samples in VCFs
    PROGRESS("Collecting samples")
    SAMPLE=[]
    vcf_table={}
    PROGRESS("There are %s Samples to be merged... "%str(len(vcfs)))
    PROGRESS('Begin ...  %s'%str(datetime.now().strftime('%H:%M:%S')))
    SAMPLES=[]
    for VCF in vcfs:
        Spl = vcf.Reader(filename=VCF).samples[0]
        if Spl in os.listdir(tmp):
            ##df=pd.read_csv(tmp+Spl, sep='\t', low_memory=False)
            #PROGRESS('****Sample %s has been formatted already. Moving on...'%Spl )
            SAMPLE.append(Spl)
            continue
        with open(VCF, 'r') as f:
            lines = [l for l in f if not l.startswith('##')] 
            Head = [l for l in f if l.startswith('##')]
               
        Table = pd.read_table( io.BytesIO(str.join(os.linesep, lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}).rename(columns={'#CHROM': 'CHROM'})           
        n=Table.shape[0] ###len() 
        PROGRESS('**VCF Opened ... %s'%str(datetime.now().strftime('%H:%M:%S')))
        Alt = ["."]*n
        Qual= ["."]*n
        Fil = ["."]*n
        Fmt = ["GT:GB"]*n 
        PROGRESS('Qual and Alt created ... %s'%str(datetime.now().strftime('%H:%M:%S')))
        Pos = list(Table["POS"])
        CH  = list(Table["#CHROM"]) 
        Ref = list(Table["REF"])
        Id  = list(Table["ID"])
        PROGRESS('Ref and ID... done... %s'%str(datetime.now().strftime('%H:%M:%S')))
        Info = [';'.join(['PERIOD='+str(record.INFO['PERIOD']), 'START='+str(record.INFO['START']), 'END='+str(record.INFO['END'])]) for record in vcf.Reader(filename=VCF)]  #**
        Unit = [x.split(';')[0] for x in Info]
        PROGRESS('INFO and Unit size done... data in frame..%s'%str(datetime.now().strftime('%H:%M:%S')))
        geno = [Getgb(record) for record in vcf.Reader(filename=VCF) ]  #**
        Spl = list(Table.columns)[-1]
        #geno = [':'.join( x.split(':')[:3] ) for x in list(Table[Spl])] #Try this for speed improve
        End = [record.INFO['END'] for record in vcf.Reader(filename=VCF)] #**
                   
        df = pd.DataFrame({"#CHROM":CH,"POS":Pos,"ID":Id,"REF":Ref ,"ALT":Alt,"QUAL":Qual,"FILTER":Fil, "INFO":Info,"FORMAT":Fmt, Spl:geno, 'U.size':Unit , 'END':End})
        SAMPLE.append(Spl)
        df[attrib+[Spl,"U.size","END"]].to_csv(tmp+Spl,sep='\t', index=None)
 
    # Output
    ##Merging
    PROGRESS(".\nMerging samples\n%s"%str(datetime.now().strftime('%H:%M:%S')))
    Y=[]
    for C in CHR: 
        if C+"_Merged.tsv" in os.listdir(tmp):
            Y.append(pd.read_csv(tmp+C+"_Merged.tsv", sep='\t'))
            PROGRESS(".\nMerged file exist %s"%str(datetime.now().strftime('%H:%M:%S')))
            continue
        Table={}
        for Spl in SAMPLE:
            df=pd.read_csv(tmp+Spl, sep='\t', low_memory=False)
            df['#CHROM'] = df['#CHROM'].astype(str)
            Table[Spl] = df.loc[df['#CHROM']==C]
        M = Mergesamples(Table)
        M['REF']=M['combined_ref']
        M['INFO']=M['combined_info']
        M =M[attrib+SAMPLE+['U.size','END']].sort_values(['#CHROM','POS'])
        PROGRESS(".\nCHR.. done %s"%str(datetime.now().strftime('%H:%M:%S')))
    ##Remove homopolymers
        if HOMO_POLY:
            M = removehomopolymers(M)
            PROGRESS("\tHomopolymers removed in data.. %s"%str(datetime.now().strftime('%H:%M:%S')))
        M.to_csv(tmp+C+"_Merged.tsv",sep='\t', index=None)
        PROGRESS(".\n..........Merged chr %s"%str(C))
        Y.append(M)
    Merged_f = pd.concat(Y) 
    
    #Format - Select columns
    
    PROGRESS("Format and sort...%s"%str(datetime.now().strftime('%H:%M:%S')))
    Merged_f = Merged_f[attrib+SAMPLE]
    
    ## Write VCF
    PROGRESS("Saving file... %s"%str(datetime.now().strftime('%H:%M:%S')))
    Merged_f.to_csv("table.tab", sep='\t', index=None)
    command = "grep '^##' "+vcfs[0]
    vcfheader = subprocess.check_output(command, shell=True)
    f=open('tmp','w')
    f.write(vcfheader.decode('utf-8'))
    f.close()
    command = "cat tmp table.tab >"+OUT+OUTFILE 
    MG = subprocess.check_output(command, shell=True)
    command = "rm table.tab"
    MG = subprocess.check_output(command, shell=True)
    #compress and index vcf 
    PROGRESS("Indexing")
    command = "bgzip -c "+ OUT+OUTFILE  +" > "+OUT+OUTFILE +'.gz'
    MG = subprocess.check_output(command, shell=True)
    command = "tabix -p vcf "+ OUT+OUTFILE  +'.gz '
    MG = subprocess.check_output(command, shell=True)
'''
 END                 
            
        vcf_reader = vcf.Reader(filename=VCF)
        PROGRESS('**VCF Opened ... %s'%str(datetime.now().strftime('%H:%M:%S')))
**        geno = [Getgb(record) for record in vcf.Reader(filename=VCF) ]
        Pos = [record.INFO['START'] for record in vcf.Reader(filename=VCF)] #
        PROGRESS('Start and Genotypes done...  %s'%str(datetime.now().strftime('%H:%M:%S')))
        CH = [record.CHROM for record in vcf.Reader(filename=VCF)]#
        PROGRESS(' CH done    %s'%str(datetime.now().strftime('%H:%M:%S')))
        Ref = [record.REF for record in vcf.Reader(filename=VCF)]#
**        End = [record.INFO['END'] for record in vcf.Reader(filename=VCF)]
        Alt = ["."]*len(CH) #
        Qual= ["."]*len(CH) #
        Fil = ["."]*len(CH) #
        Fmt = ["GT:GB"]*len(CH)#
        PROGRESS('Ref and ... done... %s'%str(datetime.now().strftime('%H:%M:%S')))
        Id = [record.ID for record in vcf.Reader(filename=VCF)] #
        PROGRESS(' ID done  %s'%str(datetime.now().strftime('%H:%M:%S')))
**        Info = [';'.join(['PERIOD='+str(record.INFO['PERIOD']), 'START='+str(record.INFO['START']), 'END='+str(record.INFO['END'])]) for record in vcf.Reader(filename=VCF)] ###
        PROGRESS( 'Info done... %s'%str(datetime.now().strftime('%H:%M:%S')))
**        Unit = [str(record.INFO['PERIOD']) for record in vcf.Reader(filename=VCF)] #
        PROGRESS('INFO and Unit size done... data in frame.. %s'%str(datetime.now().strftime('%H:%M:%S')))
    
    
 
'''
    
    
