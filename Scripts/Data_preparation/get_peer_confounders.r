#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(peer)
library(impute)


##USAGE ./path/to/get_peer_confounders.r expression_file_name
##
##	It is better to create a DIR for each tissue type & call this code from that DIR
##	This operation makes sure all output files are created in the RUNDIR
##DEFAULT is WHoleBlood

# Get argument
if (length(args)==0) {
  expr_st<- read.table('/storage/szfeupe/data/Expression_by_Tissue/Whole-Blood.rpkm', sep=" ", header=FALSE)  
} else if (length(args)==1) { 
  expr_st <- read.table(args[1], sep="\t", header=TRUE)	
}

# Set input
dim(expr_st)
#expr_st$V1 = NULL
head(colnames(expr_st))
rownames(expr_st) <-expr_st$Name
expr_st$Name = NULL
dim(expr_st)

# Remove rows of median 0 values
expr_st$rmeds <- apply(expr_st, 1, median)
expr_set<-expr_st[!(expr_st$rmeds== 0),]
expr_set$rmeds = NULL
paste("Median removal step corrected")
dim(expr_set)

# Observing the data and remove the only Asian Genome from peer analysis
ids = colnames(expr_set)
write.table(ids, file="Ids.txt")
colnames(expr_set)<-ids

dim(expr_set)

# Saving the table for EXPRANNOT
write.table(expr_set, file="Clean_expression.tsv", sep="\t")

# Set K
k = ncol(expr_set)*0.25        #k=15	
if(k > 100) {
      k <- 100 }

model = PEER()

# Set number of “hidden factors” searched for
PEER_setNk(model,k)
message("Number of factors ... ",k)

# PEER ask NxG matrix, where N=samples and G=genes
PEER_setPhenoMean(model, t(as.matrix(expr_set)))
PEER_setAdd_mean(model, TRUE)
message("Number of factors ... ",k)
head(colnames(expr_set))

# These can be changed to the tolerance you desire
PEER_setTolerance(model, .001)
PEER_setVarTolerance(model, 0.0005)

# Set limit number of iteration... most cases don’t go over 200 iterations        
PEER_setNmax_iterations(model, 250)

# Run peer on the model
PEER_update(model)

# Observing the results
factors = PEER_getX(model)
dim(factors)
weights = PEER_getW(model)
dim(weights)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
dim(residuals)
plot(precision)
rownames(factors)<- colnames(expr_set)

# Write output to files
write.table(residuals, file="peerResiduals.tsv", sep="\t")
write.table(factors, file="peerFactors.tsv", sep="\t")
write.table(weights, file="peerWeights.tsv", sep="\t")

q()
