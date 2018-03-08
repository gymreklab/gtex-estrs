#!/usr/bin/env Rscript

library(preprocessCore)

Normalize <- function(m){ apply(m,2,function(x){(x - mean(x))/sd(x)} ) } #Using normalize.quantiles.use.target instead
##	Function that generate the Hat Matrix
HatMatrix <- function(c){ return (c %*% solve(t(c) %*% c, tol=1e-23) %*% t(c))}
##	Function that generates the residuals matrix
GenerateResidual <- function(h,y){ return((diag(nrow(h)) - h) %*% y)}


# Building covariate matrix

##	 Open covariate files & set up variables
# pcas <- read.table("/storage/szfeupe/Runs/GTEx_estr/gtex.pca", sep=" ", header=FALSE)**
pcas <- read.table("/storage/szfeupe/Runs/GTEx_estr/gtex650.pca", sep=" ", header=FALSE)
rownames(pcas) <- pcas$V1
pcas$V1 = NULL
peers <- read.table("peerFactors.tsv", sep="\t", header=TRUE)
rownames(peers) <- gsub("[[:punct:]]", "-", rownames(peers))

##	 Concat PCs and PEERs
Covar <- merge(peers, pcas, by="row.names")
rownames(Covar) <- Covar$Row.names
Covar$Row.names = NULL
dim(Covar)


#	 Load Expression data
Expr <- read.table("Clean_expression.tsv", sep="\t", header=TRUE)
colnames(Expr) <- gsub("[[:punct:]]" , "-", colnames(Expr))
dim(Expr)
Yo=Expr[,rownames(Covar)]
paste("Covariates [check] ... Expression [check]")

##	Function to normalize expression: For clean_Expr only (We don't use this)


##	Transpose expression matrix to get samples in rows
Y <- data.matrix(t(Yo))
dim(Y)

##	Test
hist(Y[,"ENSG00000237683.5"])
hist(Y[,"ENSG00000188976.6"])


#Fitting expression into normal dist mean=0, sd=1
Y1=normalize.quantiles.use.target(as.matrix(Y), rnorm(nrow(Y)))
dimnames(Y1) = dimnames(Y)
dim(Y1[,1])
hist(Y1[,"ENSG00000237683.5"])
hist(Y1[,"ENSG00000188976.6"])
#write.table(Y1, file="TT.Clean_expression.csv", sep=",")     ##Not needed

## 	Transform covariate to matrix and call it C
C <- data.matrix(Covar)
dim(C)
paste(" Covariates Formatting [check]")


##	CORRECTION
#	Let's get the Hat matrix
H <- HatMatrix(C)
dim(H)

## 	Let's get the residual expression matrix Y.prime 
Y.prime <- GenerateResidual(H,Y1)
dim(Y.prime)
#head(rownames(Y.prime))

##	Record it
write.table(Y.prime, file="Corr_Expr.csv", sep=",")

paste( " Corrected Expression (Residuals) [check].... Adding pop...")

# We can ADD population covariates, here: age and sex but there are many missing data
#
#
q()
