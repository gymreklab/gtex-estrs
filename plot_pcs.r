#!/usr/bin/Rscript

########USAGE
#./plot_pcs pcfile1 pcfile2
#pcfile2 is optional . The file have at least 3 columns (PC1,PC2,PCC) PCC=population label

#################This section may be commented to enter file within the code on args[1] alone or with args[2]
#Capture one or 2 file with PCAs
args = commandArgs(trailingOnly=TRUE)
#The default if no argument is entered
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

paste("File1 is", args[1])
paste("File1 is", args[2])

#names(df)=c("IDs","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PCC")
names(df)
#head(df)
 
#plot both with different colors (GTEx in shades of gray)
plot(df$PC1, df$PC2,xlab="PC1",ylab="PC2", col = ifelse(df$PCC=='AMR','coral1', ifelse(df$PCC=='AFR','darkseagreen1',ifelse(df$PCC=='EUR', 'yellow', ifelse(df$PCC=='EAS','khaki1', ifelse(df$PCC=='SAS','lightblue', ifelse(df$PCC=='AfricanAmerican','gray87', ifelse(df$PCC=='Asian','gray60','gray40'))))))), pch=20)
legend("topleft",title="Race",legend=c('AFR','African American','AMR','EAS', 'SAS','Asian','EUR','European'),pch=20, col=c('darkseagreen1','gray87','coral','khaki1','lightblue','gray60','yellow','gray40'))
title(main="PCA 1000G & GTEx")

#Plot 1000G only
#plot(df$PC1, df$PC2,xlab="PC1",ylab="PC2", col = ifelse(df$PCC=='AMR','coral1', ifelse(df$PCC=='AFR','darkseagreen1',ifelse(df$PCC=='EUR', 'dimgray', ifelse(df$PCC=='EAS','khaki1', ifelse(df$PCC=='SAS','lightblue', ifelse(dfs$PCC=='AfricanAmerican','green', ifelse(dfs$PCC=='Asian','blue','gray'))))))), pch=20)
#title(main="PCA 1000G")

#plot both with 1000G in gray scale
#plot(df$PC1, df$PC2,xlab="PC1",ylab="PC2", col = ifelse(df$PCC=='AMR','gray87', ifelse(df$PCC=='AFR','gray72',ifelse(df$PCC=='EUR', 'dimgray', ifelse(df$PCC=='EAS','gray58', ifelse(df$PCC=='SAS','gray16', ifelse(dfs$PCC=='AfricanAmerican','green', ifelse(dfs$PCC=='Asian','blue','gray'))))))), pch=20)
#legend("topleft",title="Race",legend=c('AFR','African American','AMR','EAS', 'SAS','Asian','EUR','White'),pch=20, col=c('gray72','green','gray87','gray58','gray16','red','dimgray', 'blue') )
title(main="PCA 1000G & GTEx")

#Plot GTEx pcs only
#plot(dfs$PC1, dfs$PC2,xlab="PC1",ylab="PC2", col = ifelse(dfs$PCC=='AfricanAmerican','green', ifelse(dfs$PCC=='Asian','blue','blue')), pch=20)
#legend("topleft",title="Race",legend=c('Asian','African American','European'),pch=20, col=c('red','green','blue') )
#title(main="GTEx only")

paste("... Done")
quit()


