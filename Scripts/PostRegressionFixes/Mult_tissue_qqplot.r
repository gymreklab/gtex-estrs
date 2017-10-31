#!/usr/bin/env Rscript

## Get the data
args = commandArgs(trailingOnly=TRUE)
#Default if no argument is entered
if (length(args)==0) {
  args[1]="Lin_Reg_Out"
  args[2]="Lin_Reg_Out_Permute"
  		     }
 
#wb="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/WholeBlood/Lin_Reg_OutFin.txt"
cf="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/Cells-Transformedfibroblasts/Lin_Reg_OutFin.txt"
ms="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/Muscle-Skeletal/Lin_Reg_OutFin.txt"
at="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/Artery-Tibial/Lin_Reg_OutFin.txt"
as="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/Adipose-Subcutaneous/Lin_Reg_OutFin.txt"
lg="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/Lung/Lin_Reg_OutFin.txt"
em="/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/Esophagus-Mucosa/Lin_Reg_OutFin.txt"

linout = read.table(args[1], sep="\t", header=TRUE)
permuted = read.table(args[2], sep="\t", header=TRUE)

linout1 = read.table(cf, sep="\t", header=TRUE)
linout2 = read.table(ms, sep="\t", header=TRUE)
linout3 = read.table(at, sep="\t", header=TRUE)
linout4 = read.table(as, sep="\t", header=TRUE)
linout5 = read.table(lg, sep="\t", header=TRUE)
linout6 = read.table(em, sep="\t", header=TRUE)
# PVALUES
logpval = -log(linout$p.wald, base=10)
logpva2 = -log(linout1$p.wald, base=10)
logpva3 = -log(linout2$p.wald, base=10)
logpva4 = -log(linout3$p.wald, base=10)
logpva5 = -log(linout4$p.wald, base=10)
logpva6 = -log(linout5$p.wald, base=10)
logpva7 = -log(linout6$p.wald, base=10)
#Uniform disttribution
Unifd = runif(length(logpval) , min=0, max=1)
Unif_dist = -log(Unifd, base=10)
head(Unif_dist)
#QQ plot logpval vs inuf.dist
png(filename="QQplot2.png")
##qqplot(logpval, Unif_dist, col='blue',  xlab = "Observed Qunantiles", ylab = "Theoretical Quantile")
##qqplot(Unif_dist, logpval, col='blue')
A <- qqplot(Unif_dist, logpval, plot.it=FALSE )
B <- qqplot(Unif_dist, logpva2, plot.it=FALSE )
C <- qqplot(Unif_dist, logpva3, plot.it=FALSE )
D <- qqplot(Unif_dist, logpva4, plot.it=FALSE )
E <- qqplot(Unif_dist, logpva5, plot.it=FALSE )
F <- qqplot(Unif_dist, logpva6, plot.it=FALSE )
qqplot(Unif_dist, logpva7, col='orange',  ylab = "Observed Qunantiles", xlab = "Theoretical Quantile")
paste("Observed Pvalues on  Uniform distribution")
##Permuted
logper = -log(permuted$p.wald, base=10)
#QQ plot logpermuted vs logobserved
Permut_pts <- qqplot(Unif_dist, logper, plot.it=FALSE)
points(Permut_pts, col="grey")
points(A, col="red")
points(B, col="black")
points(C, col="yellow")
points(D, col="pink")
points(E, col="brown")
points(F, col="cyan")

legend("topleft",legend=c("Whole blood","Fibroblast","Muscle Skeletal","Arterial Tibial","Adipose Subcutaneous","Lung","Esophagous mucosa", "Permuted control"), col=c("red","black","yellow","pink", "brown","cyan","orange","grey"), pch=c(1,3))

paste("Permutted Pvalues ovelayed on Observed Pvalues")

#Trace Diagonal for better separation
lines(x = c(0,8), y = c(0,8))
dev.off()

## save data into a file
#write.table(linout, file="Adjusted_Lin_Reg_Out", sep="\t")


q() 

