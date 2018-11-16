#!/usr/bin/env python
import pandas as pd
import numpy as np
import scipy.stats as ss
import argparse
import sys

"""
USAGE: ./merge_ld_causal.py <LDFILE> <OUTFILE><CATALOG>
CATALOG is optinal

"""

GWASCAT='/storage/szfeupe/Runs/650GTEx_estr/gwas/gwas_catalog_loci_all.tab' 
CAUSAL = '/storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Merged_Best_causality.Table'
GWASLD = '/storage/szfeupe/Runs/650GTEx_estr/gwas/gwas_catalog_loci_all.tab'

try:
    LD01 = sys.argv[1]
    OUTFILE = sys.argv[2]
except:
    LD01 = '/storage/szfeupe/Runs/650GTEx_estr/gwas/str_snp_ld_01_up.tab'
    OUTFILE = 'GWAS_01.tab'
    print('Using this file: ', LD01, ' and output in ', OUT)
    sys.stderr.write(__doc__)

try: GWASCAT = sys.argv[3]
except: print("Using catalog ", GWASCAT)

catalog = pd.read_csv(GWASCAT, sep='\t')
catalog['snp_locus'] = catalog.apply(lambda x: str(x['chrom'])+':'+str(int(x['snp.start'])) , 1)
#catalog =catalog[['chrom','snp.start','rsid','locusname','pval', 'snp_locus']].copy()
#
ld_data = pd.read_csv(LD01, sep='\s+').drop_duplicates().dropna()
ld_data['pval'] = ld_data['pval'].astype(float)
ld_data.columns = ['str_locus', 'snp_locus','allele','freq_het','MAF','KL','r2','pval']
gwas_loci = ld_data.loc[(ld_data['pval'].astype(float)>0)&(ld_data['pval'].astype(float)<0.05)]

#
trait_map = pd.read_csv(GWASLD, sep='\t')
trait_map['snp_locus'] = trait_map.apply(lambda x: x['chrom']+':'+str(x['snp.start']) , 1)
#
estrs = pd.read_csv(CAUSAL, sep='\t').dropna()
estrs['str_locus'] = estrs.apply(lambda x: x['chrom'][3:]+':'+str(int(x['best.str.start'])) , 1)
estrs = estrs.sort_values('str_locus')
estrs['best.score'] = estrs['best.score'].apply(lambda x: 100*x , 1)
#
#Merge to get GWAS1
gwas_str = gwas_loci.merge(estrs, on =['str_locus'], how='inner')
gwas_str.to_csv(OUTFILE, sep='\t', index=None)
print list(ld_data.columns),gwas_str.shape, '\n',list(gwas_str.columns),'\n', list(catalog.columns)
#gwas_str.merge(catalog.drop_duplicates()[['rsid','locusname','snp_locus']], on=['snp_locus'], how='inner').to_csv(OUTFILE, sep='\t', index=None)


