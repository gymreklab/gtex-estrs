import argparse
import pandas as pd
import numpy as np
import scipy.stats as ss
import sys
"""
USAGE: getld.py <SNP_LOCI_FILE> <STR_LOCI_FILE>
SNP_LOCI_FILE formatted as csv or tab delomited with 3 mandatory columns chrom, snp_pos, gene or gene.name
"""
STRGT = '/storage/szfeupe/Runs/650GTEx_estr/Genotypes/NormalizedGenotypes.table'
SNPGT = '/storage/szfeupe/Runs/650GTEx_estr/SNP_Analysis/'        #%chr??.tab
DISTFROMrsid = 10000 #10kb


args = parser.parse_args()
gwas = pd.read_csv(args[0], sep='\t',low_memory=False ).sort_values("chrom")
print (gwas.shape[0],' gwas loci entered')

loci_str = pd.read_csv(STRGT, sep='\t')

for chrom in list(gwas['chrom']):
    pd.read_csv('%s')

for i in range(expr_annot.shape[0]):
        gene = expr_annot.index.values[i]
        PROGRESS(" Getting data for %s"%gene)
        start = expr_annot["gene.start"].values[i]
        end = expr_annot["gene.stop"].values[i]
        
        PROGRESS("%s STRs tested \n"%str(cis_strs.shape[0]))
       












# genotypes
strgt = pd.read_csv(STRGT, sep='\t')

#Overlap by gene
if 'gene' in gwas.index :
    trait = gwas.loc[gwas['gene'].isin(list(caus_estrs['gene']))].merge(caus_estrs, on=['chrom','gene'])
elif 'gene.name' in gwas.index : 
    trait = gwas.loc[gwas['gene.name'].isin(list(caus_estrs['gene.name']))].merge(caus_estrs, on=['chrom','gene.name'])
else:
    print ('Error in SNP file... No gene column to make connection...')
    sys.exit(1)
    
trait['locus'] = trait.apply(lambda x: x['chrom']+':'+str(x['str.start']),1)

print('Overlaps on gene with our study: ',trait.shape)
result={}
result['locus']=["gene","snp_start", "str.start", "Tissue", "num.samples","r_square","ld.pvalue"]
# Calculate LD
for chrom in set(trait['chrom']):
    print('Starting with chrom ',chrom)
    data = trait.loc[trait['chrom']==chrom]
    if data.shape[0]==0: print('No data for %s...'%chrom); continue
    snpgt = pd.read_csv('%s%s.tab'%(SNPGT,chrom), sep='\t',low_memory=False)
    gtstr = strgt.loc[strgt['chrom']==chrom]
    print('There is data in %s'%chrom, data.shape)
    for index, row in data.iterrows():
        snp_start=int(row['snp.start'])
        str_start=int(row['str.start'])
        tissue = row['best.tissue']
        locus = row['locus']
        gene = row['gene']
        pheno = row['traits']
        genename=row['gene.name']
        Xstr = gtstr.loc[gtstr['start']==str_start]
        if Xstr.shape[0]==0: print(str_start,' passed'); continue  # No str... This should not happend
        Xsnp = snpgt.loc[snpgt['start']==snp_start]
        if Xsnp.shape[0]==0: print(snp_start,' passed'); continue  # No snp
        samples = list(pd.read_csv('%s%s/Corr_Expr.csv'%(EXPRESSION,tissue)).index)
        s = pd.concat([Xstr[samples].replace('None', np.nan),Xsnp[samples].replace('None', np.nan)])
        s = s.astype(float).dropna(axis=1)
        ld_value = ss.pearsonr(s.loc[s.index[0],:], s.loc[s.index[1],:])
        result[locus]=[gene,snp_start, str_start, tissue, len(samples),ld_value[0],ld_value[1]]
        
df=pd.DataFrame.from_dict(result)
df1 = df.transpose()
df1.columns = result['locus']
df1.to_csv('LD_FILE.tsv',sep='\t')
print 'ld file created as LD_FILE.tsv' 