#!/usr/bin/env Rscript

# Peer Factors = peerFactors.tsv
# IDs for Peer Factors in order = head -1 Clean_expression.tsv
# PCs factors =  /home/szfeupe/projects/GTEX_estrs/data/gtex_pca_only
# Age, Gender, Race = Phen_attributes
# awk 'FNR==NR {a[$1]=$0; next}; $1 in a {print a[$1],$2,$3}' gtex.pca phen_attributes

library(preprocessCore)


# Building covariate matrix

##	 Open covariate files & set up variables

pcas <- read.table("/storage/szfeupe/Runs/GTEx_estr/gtex.pca", sep=" ", header=FALSE)
rownames(pcas) <- pcas$V1
pcas$V1 = NULL
peers <- read.table("peerFactors.tsv", sep="\t", header=TRUE)
rownames(peers) <- gsub("[[:punct:]]", "-", rownames(peers))


##	 Concat them
Covar <- merge(peers, pcas, by="row.names")
write.table(Covar, file="Merged.tsv", sep="\t")
#head(rownames(Covar))
rownames(Covar) <- Covar$Row.names
Covar$Row.names = NULL
dim(Covar)


##	 Load Expression data
Expr <- read.table("Clean_expression.tsv", sep="\t", header=TRUE)
colnames(Expr) <- gsub("[[:punct:]]" , "-", colnames(Expr))
dim(Expr)
#head(Expr[,1:5])
Yo=Expr[,rownames(Covar)]
#head(Yo[,1:5])
paste("Covariates [check] ... Expression [check]")


##	Function to normalize expression: For clean_Expr only (Not used here)
Normalize <- function(m){ apply(m,2,function(x){(x - mean(x))/sd(x)} ) }

##	Function that generate the Hat Matrix
HatMatrix <- function(c){ return (c %*% solve(t(c) %*% c, tol=1e-23) %*% t(c))}

##	Function that generates the residuals matrix
 GenerateResidual <- function(h,y){ return((diag(nrow(h)) - h) %*% y)}

##	Transpose expression matrix to get samples in rows
Y <- data.matrix(t(Yo))		
dim(Y)


##	Test
hist(Y[,"ENSG00000237683.5"])
hist(Y[,"ENSG00000188976.6"])


#write.table(Y, file="T.Clean_expression.csv", sep=",")
Y1=normalize.quantiles.use.target(as.matrix(Y), rnorm(143))
dimnames(Y1) = dimnames(Y)
dim(Y1[,1])
hist(Y1[,"ENSG00000237683.5"])
hist(Y1[,"ENSG00000188976.6"])
#write.table(Y1, file="TT.Clean_expression.csv", sep=",")


## 	Transform covariate to matrix and call it C
C <- data.matrix(Covar)
dim(C)
paste(" Functions and Formatting [check]")


##	CORRECTION
##	Let's get the Hat matrix
H <- HatMatrix(C)
dim(H)

## 	Let's get the residual expression matrix Y.prime 
Y.prime <- GenerateResidual(H,Y)
dim(Y.prime)
#head(rownames(Y.prime))

##	Record it
write.table(Y.prime, file="Corr_Expr.csv", sep=",")

paste( " Corrected Expression (Residuals) [check].... Adding pop...")
 
# We can ADD population covariates, here: age and sex but there are many missing data

##	New Covariate matrix
P = read.table('/storage/szfeupe/Runs/GTEx_estr/phen_attributes', sep='\t', header=TRUE, na.strings=c("Donor","(OPO)"))
						 #[SUBJID-GENDER-AGE-RACE-ETHNCTY]
rownames(P) <- P$SUBJID
P$SUBJID <- NULL

####################out with race and ethnicity. Was added in PCA table
P$ETHNCTY = NULL
P$RACE=NULL
####################
pheno=P
paste(" ...  POP COVARIATES") #print(pheno, topn=10)
dim(pheno)
dim(Covar)
covariate = merge(Covar, pheno, by="row.names")

dim(covariate)
head(covariate[,25:30])
paste(" ...  Pop Residuals")
##

covari<- data.matrix(covariate)
Hpop <- HatMatrix(covari)
Y.prime.pop <-GenerateResidual(Hpop, Y)
write.table(Y.prime.pop, file="Corr_with_pop_Expr.csv", sep=",")


q()
