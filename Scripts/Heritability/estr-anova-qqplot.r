#!/usr/bin/env Rscript

## Get the data
args = commandArgs(trailingOnly=TRUE)
#Default if no argument is entered
if (length(args)==0) {
  args[1]="Lin_Reg_Out"
  args[2]="Lin_Reg_Out_Permute"
  		     }

#Get estr anova table

linout = read.table(args[1], sep="\t", header=TRUE)
#permuted = read.table(args[2], sep="\t", header=TRUE)


# PVALUES
logpval = -log(linout$anova_pval, base=10)

#Uniform disttribution
Unifd = runif(length(logpval) , min=0, max=1)
Unif_dist = -log(Unifd, base=10)
head(Unif_dist)

#QQ plot logpval vs inuf.dist
png(filename="QQplot0.png")


#qqplot(Unif_dist, logpval, col=factor(linout$strprefer),  ylab = "Observed P values (-log10)", xlab = "Theoretical P values (-log10)")

qqplot(Unif_dist, logpval, col='blue',  ylab = "Observed P values (-log10)", xlab = "Theoretical P values (-log10)")


paste("Observed Pvalues on  Uniform distribution")

##Permuted
logper = -log(Unifd, base=10)

#QQ plot logpermuted vs logobserved
Permut_pts <- qqplot(Unif_dist, logper, plot.it=FALSE)
#print(paste(Permut_pts))
points(Permut_pts, col="red")

legend("topright",legend=c("Observed P-values", "Negative control"), col=c("blue","red"), pch=c(1,3))

paste("Permutted Pvalues ovelayed on Observed Pvalues")

#Trace Diagonal for better separation
lines(x = c(0,5), y = c(0,5))
dev.off()

q() 

