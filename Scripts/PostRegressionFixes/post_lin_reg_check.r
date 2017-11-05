#!/usr/bin/env Rscript

## Get the data
args = commandArgs(trailingOnly=TRUE)
#Default if no argument is entered
if (length(args)==0) {
  args[1]="Lin_Reg_Out"
  args[2]="Lin_Reg_Out_Permute"
  		     }

linout = read.table(args[1], sep="\t", header=TRUE)
permuted = read.table(args[2], sep="\t", header=TRUE)


# PVALUES
logpval = -log(linout$p.wald, base=10)

#Uniform disttribution
Unifd = runif(length(logpval) , min=0, max=1)
Unif_dist = -log(Unifd, base=10)
head(Unif_dist)

#QQ plot logpval vs inuf.dist
png(filename="QQplot0.png")
##qqplot(logpval, Unif_dist, col='blue',  xlab = "Observed Qunantiles", ylab = "Theoretical Quantile")
qqplot(Unif_dist, logpval, col='blue',  ylab = "Observed Qunantiles", xlab = "Theoretical Quantile")

paste("Observed Pvalues on  Uniform distribution")

##Permuted
logper = -log(permuted$p.wald, base=10)

#QQ plot logpermuted vs logobserved
Permut_pts <- qqplot(Unif_dist, logper, plot.it=FALSE)
points(Permut_pts, col="grey")
legend("topright",legend=c("Original data", "Permuted"), col=c("blue","grey"), pch=c(1,3))

paste("Permutted Pvalues ovelayed on Observed Pvalues")

#Trace Diagonal for better separation
lines(x = c(0,8), y = c(0,8))
dev.off()

## save data into a file
#write.table(linout, file="Adjusted_Lin_Reg_Out", sep="\t")


q() 
