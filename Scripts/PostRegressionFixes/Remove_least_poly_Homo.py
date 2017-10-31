tissu = ['WholeBlood','Cells-Transformedfibroblasts','Muscle-Skeletal','Lung','Adipose-Subcutaneous', 'Artery-Tibial','Esophagus-Mucosa']

##This code will remove homopolymer and least polymorphic sitesfrom linear regression analysis output
import pandas as pd
import numpy as np
import sys

args = parser.parse_args()
Tissu = args.tissue

##It will remove locus that fall in homopolymer or least polymorphic sites from linear regression analysis output

def filterregression(Tissu):
    #
    print(Tissu)
    chrs = [i for i in range(1,23)]
    chrs.append('X') ; chrs.append('Y')

    #regression
    Input='/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/'+Tissu
    linreg = pd.read_csv(Input+'/Lin_Reg_Out', sep='\t')
    #linreg = pd.read_csv(Input+'/Lin_Reg_Out_perm', sep='\t')

    #Homopolymers
    strfile= "/storage/szfeupe/Runs/GTEx_estr/Homopolymer_sites.tsv"
    L_pol = pd.read_csv(strfile, sep='\t')

    index=['0']
    Out = pd.DataFrame(index=index, columns=linreg.columns)

    #LR Cleanup
    for i in chrs:
    #Open sites to be removed
        Lpolym = '/storage/szfeupe/Runs/GTEx_estr/Heterozygosity/Heterozygosity_VCF.chr'+str(i)
        Reg_to_remove=pd.read_csv(Lpolym, sep='\t', header=None)
        Reg_to_remove.columns=['CHROM','START','END','NumAlleleLenghts','NumAlleles','Het.length','Het.genotype']
        Homopoly = L_pol.loc[L_pol['chrom']=='chr'+str(i)]

    #   Remove Homoploymers
        LIN = linreg.loc[linreg['chrom']=='chr'+str(i)]
        LIN0 = LIN.loc[LIN['str.start'].isin(list(Homopoly['str.start']))==False]

    #   Remove least polymorphic
        Reg=Reg_to_remove.loc[Reg_to_remove['Het.length']<=0.3]
        FIN = LIN0.loc[LIN0['str.start'].isin(list(Reg['START']))==False]

    #   Sum it up and outout   
        frames = [Out, FIN]
        Out = pd.concat(frames)

        print(LIN.shape, LIN0.shape, FIN.shape, 'End')

    Out1=Out.drop('0')
    Out1=Out1.drop('Unnamed: 0', 1)
    Out1.to_csv(Input+'/Lin_Reg_OutFin.txt', sep='\t', index=False)

#for T in tissu:
#    filterregression(T)
filterregression(Tissu)