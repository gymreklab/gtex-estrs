##This code creates files from the VCF file containing least polymorphic sites and homopolymers
##This file will be used to filter sites after linear regression

import pandas as pd
import numpy as np
import vcf

chrs = range(1,23,1)
chrs.append('X') ; chrs.append('Y')


Out=open('/storage/szfeupe/Runs/GTEx_estr/Homopolymer_sites.tsv', 'w')
Out.write('\t'.join(['chrom','str.start','unit.length', 'str.end'])+'\n')

#VCF filtering 
for i in chrs:
#   Collecting least polymorphic sites
    strfile= "/storage/szfeupe/Runs/GTEx_estr/STRs_No_Missing_genotypes/STR_filter.chr"+str(i)+".vcf.gz"
    STRs = vcf.Reader(filename=strfile)
    print "Chr",i,' ...'
    for record in STRs:
        if record.INFO['PERIOD'] ==1:
            Out.write("\t".join([str(record.CHROM), str(record.POS), str(record.INFO['PERIOD']), str(record.INFO['END'])])+'\n')
    print "Done with chr",i,' ...'
Out.close()

print "done..."