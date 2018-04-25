#!/usr/bin/Rscript


library(qvalue)
p <- scan("pvalues.txt")
qobj <- qvalue(p, fdr.level = 0.1, pfdr=TRUE) #, pi0.method="bootstrap" )
write.qvalue(qobj,file="qvalues.txt")

#paste (qobj) 

#summary(qobj)
hist(qobj)
plot(qobj)

print('the end')

