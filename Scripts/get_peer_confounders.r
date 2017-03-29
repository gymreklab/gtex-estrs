#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(peer)
library(impute)

##USAGE ./path/to/get_peer_confounders.r expression_file_name
##
##	It is better to create a DIR for each tissue type & call this code from that DIR
##	This operation makes sure all output files are created in the RUNDIR
##DEFAULT is WHoleBlood

#Get argument
if (length(args)==0) 
{  expr_set<- read.table('/storage/szfeupe/data/Expression_by_Tissue/Whole-Blood_rpkm', sep=" ", header=FALSE)  } 
else if (length(args)==1) 
{ expr_set <- read.table(args[1], sep=" ", header=FALSE)	}

# Set input
dim(expr_set)
expr_set$V1 = NULL

#Remove rows of mean 0 values
expr_set <- expr_s[rowSums(expr_s[, -1] > 0) != 0, ]

#saving the table for EXPRANNOT
write.table(expr_set, file="Clean_expression.tsv", sep="\t")

#Observing the data
ids = expr_set$V1
write.table(ids, file="Clean_Exp_table/ids")
#rownames(expr_set) <-expr_set$V1
expr_set$GeneID = NULL
dim(expr_set)

#set K
k = ncol(expr_set) * .25
model = PEER()

# set number of “hidden factors” searched for
PEER_setNk(model,k)

#PEER ask NxG matrix, where N=samples and G=genes
PEER_setPhenoMean(model, t(as.matrix(expr_set)))
PEER_setAdd_mean(model, TRUE)

# These can be changed to the tolerance you desire
PEER_setTolerance(model, .001)
PEER_setVarTolerance(model, 0.0005)

# Set limit number of iteration... most cases don’t go over 200 iterations        
PEER_setNmax_iterations(model, 500)

#Run peer on the model
PEER_update(model)

#Observing the results
factors = PEER_getX(model)
dim(factors)
weights = PEER_getW(model)
dim(weights)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
dim(residuals)
plot(precision)

#Write output to files
write.table(residuals, file="peerResiduals.tsv", sep="\t")
write.table(factors, file="peerFactors.tsv", sep="\t")
write.table(weights, file="peerWeights.tsv", sep="\t")
write.table(precision, file="peerPrecision.tsv", sep="\t")
q()
