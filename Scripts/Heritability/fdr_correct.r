#!/usr/bin/Rscript

library(qvalue)
mypvalues <- read.table("pvalues.txt")
qobj <- qvalue(mypvalues, fdr.level = 0.05, pfdr=TRUE) #, pi0.method="bootstrap" )
write.qvalue(qobj,file="qvalues.txt")
# Visualize it 
hist(qobj)
plot(qobj)

print('the end')