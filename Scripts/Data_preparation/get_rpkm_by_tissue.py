import os
files = os.listdir('/home/szfeupe/projects/GTEX_eSTRs/data/RNA-Seq/')

	#######Head_wgs is the file with all IDs of sample with WGS

exp = open('head_wgs','r')
exp = exp.readlines()
exp = [x.strip('\n')for x in exp]
#exp = [x[1] for x in exp]
print len(exp), ' WGS'
	######### 'files' are the RNA-Seq per tissue. Each tissues named after itself

print '  '.join(files)
mypath='/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/ExpressionFiles/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct.gz'

	######### Remove the rpkm header with all sample IDs with rpkm values
indices = {}
header_rpkm = open('header','r').readlines()[2].split()
for n, i in enumerate(header_rpkm):
        try:
            indices[i].append(n)			#Save in a dictionary for later references
        except:
            indices[i] = [n]

print ' Length is ',len(header_rpkm)

print len(exp)
i=1+0
	########### We now screen through each tissue info file to (1) get our WGS ids of interest (2) get the rpkm
		###### 
os.system("echo 'Total number of samples: 148 for WGS \n'> statistics") #if i==1:
s=open('statistics','w')
s.write('Total number of samples: 148 for WGS\n')
for file in files:
#    file = 'Muscle-Skeletal'
    Out = open('test.txt','w')
    Out.write("")
    Out.close()

	### test.txt receives all lines from rnaseq-tissue specific file that pertain to a given wgs ID
	### For a given ID if there are ++ input, we take the most recent by col#11 (-k 11) in the Sra line

    for lines in exp:
	command = "grep '"+lines+"' ../projects/GTEX_eSTRs/data/RNA-Seq/"+file+"| sort -k 11| tail -1|sort -k 1 >> test.txt"
#	print command
    	os.system(command) 
    print command, "\nDONE*******\n",lines
#
#    command = "awk '{print $17}' test.txt | grep '^GTEX' > ids.txt"
    command = "awk '{if ($4=="+'"BI"'+"){$17=$17} else {$17=$18} print}' test.txt | awk '{print $17}' | grep '^GTEX' > ids.txt"
    os.system(command)
#
    columns=[]
#
    ff = open('ids.txt','r').readlines()
    print len(ff), ' IDs sample(s) in ', file
    hed=['GeneID']
    for lines in ff:
   	lines = lines.strip()
        #print lines
    	try:
	    c=indices[lines][0]
	    columns.append(c)
	    hed.append(lines)
    	except:
	    pass
    print 'Number of samples for ',file,' is ',len(columns),'\t',100*len(columns)/148,'% of samples'
#
    
    stat=file + '\t' + str(len(columns)) +'\t'+ str(100*len(columns)/148)+'%\n'
    s.write(file + '\t' + str(len(columns)) +'\t'+ str(100*len(columns)/148)+'%\n')
    columns =['$'+str(x) for x in columns]
#    
    command = "zcat "+mypath+ "| tail -56319 |awk '{print $1,"+','.join(columns)+" }'  > Rpkm"
#    print command
#
    os.system(command)
    command = "awk 'FNR==NR {a[$1]=$0; next}; $1 in a {print a[$1]}' Rpkm  protein_coding_IDs > /storage/szfeupe/data/Expression_by_Tissue/"+file+".rpkm"
    os.system(command)
    command = "sed -i '1 i "+ ' '.join(hed)+"' /storage/szfeupe/data/Expression_by_Tissue/"+file+".rpkm"
    print command
    os.system(command)
    os.system("rm Rpkm")

#End for
"""
tissue = [x.split('-')[2] for x in exp]
tissue = list(set(tissue))
print '\t'.join(tissue),'**',len(tissue), '**',len(exp)
for tis in tissue:
    out=open(tis, 'w')
    n=0
    for line in exp:
	if line.split('-')[2]==tis:
	    out.write(line+'\n')
	    n=n+1
    out.close()
    print tis,'\t',n
"""
print 'DONE'
