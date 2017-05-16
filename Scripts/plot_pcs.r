#!/usr/bin/Rscript

########USAGE
#./plot_pcs pcfile1 pcfile2
#pcfile2 is optional . The file have at least 3 columns (PC1,PC2,PCC) PCC=population label
#	PCs file source is /storage/mgymrek/gtex/genotypePCA//GTEx_1KG_merged.pca.evec


########This section may be commented to enter file within the code on args[1] alone or with args[2]
#Capture one or 2 files with PCAs
args = commandArgs(trailingOnly=TRUE)

#Default if no argument is entered
if (length(args)==0) {
  args[1]="/home/szfeupe/projects/GTEX_eSTRs/data/mapped_merged_1000G_only"
  args[2]="/home/szfeupe/projects/GTEX_eSTRs/data/mapped_merged_pca_GTEX_only"
  		     }

#setwd("~/projects/GTEX_eSTRs/Runs/")

#Open : This code assume both files have headers and they are similar
if (length(args)==1) {
  df <- read.table(args[1],header=TRUE)
} else if (length(args)>1) {
  df1 <- read.table(args[1], header=TRUE)
  df2 = read.table(args[2],header=TRUE)
  df <- rbind(df1, df2)
}

#paste("File1 is", args[1])	#Just checking
#paste("File2 is", args[2])

#If no header Uncomment line below
#names(df)=c("IDs","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PCC")

names(df)
#head(df)
 
#plot both with different colors (1000G in shades of gray)
plot(df$PC2, df$PC3,xlab="PC1",ylab="PC2", col = ifelse(df$PCC=='AMR','dimgray', ifelse(df$PCC=='AFR','gray58',ifelse(df$PCC=='EUR', 'gray', ifelse(df$PCC=='EAS','gray72', ifelse(df$PCC=='SAS','gray87', ifelse(df$PCC=='AfricanAmerican','yellow', ifelse(df$PCC=='Asian','blue','green'))))))), pch=20)
legend("topleft",title="Race",legend=c('AFR','AMR','EAS', 'SAS','EUR','African American','Asian','European'),pch=20, col=c('gray58','dimgray','gray72','gray87','grey','yellow','blue', 'green'))
title(main="PCA 1000G & GTEx")

#Plot 1000G only
#plot(df$PC1, df$PC2,xlab="PC1",ylab="PC2", col = ifelse(df$PCC=='AMR','coral1', ifelse(df$PCC=='AFR','darkseagreen1',ifelse(df$PCC=='EUR', 'dimgray', ifelse(df$PCC=='EAS','khaki1', ifelse(df$PCC=='SAS','lightblue', ifelse(dfs$PCC=='AfricanAmerican','green', ifelse(dfs$PCC=='Asian','blue','gray'))))))), pch=20)
#title(main="PCA 1000G")

#Plot GTEx pcs only
#plot(dfs$PC1, dfs$PC2,xlab="PC1",ylab="PC2", col = ifelse(dfs$PCC=='AfricanAmerican','green', ifelse(dfs$PCC=='Asian','blue','blue')), pch=20)
#legend("topleft",title="Race",legend=c('Asian','African American','European'),pch=20, col=c('red','green','blue') )
#title(main="GTEx only")

paste("... Done")
quit()


