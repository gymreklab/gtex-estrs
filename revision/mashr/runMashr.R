#!/bin/R

#------------------Load the data into two dataframes----------------
library('purrr')
#get files to load
files = list.files(path="workdir", pattern="*.table", full.names=TRUE)
shortFileNames = purrr::map(files, function(name) basename(tools::file_path_sans_ext(name)))

#load files individually
#lung.table <- read.table("Lung.table", header=TRUE, colClasses=c('character', 'factor', 'integer', 'numeric', 'numeric'))
dataFrameList = purrr::map(files, read.table, header=TRUE, colClasses=c('character', 'factor', 'integer', 'numeric', 'numeric'))

#append the file name to the columns named beta and beta.se
append = function (name) {
	function (df, ext) {
		names(df)[names(df) == name] = paste(name, ext, sep='_')
		return(df)
	}
}

dataFrameList = purrr::map2(dataFrameList, shortFileNames, append('beta'))
dataFrameList = purrr::map2(dataFrameList, shortFileNames, append('beta.se'))

#merge the dataframes
allData = purrr::reduce(dataFrameList, function(df1, df2) base::merge(df1, df2, by=c('gene', 'chrom', 'str.start')))

#move gene, str idetnifier info into rownames
rownames(allData) = purrr::pmap(allData[c(1,2,3)], paste, sep='_')
allData = allData[-(1:3)]

#separate betas and beta.ses
betas = allData[!grepl('se', colnames(allData))]
beta.ses = allData[grepl('se', colnames(allData))]

colnames(betas) = purrr::map(colnames(betas), substring, 6)
colnames(beta.ses) = purrr::map(colnames(beta.ses), substring, 9)

#-------------------------Use mashr------------------------
library(mashr)
#For overall coding pattern, look at this vignette:
#https://stephenslab.github.io/mashr/articles/intro_mash.html

#1)hand the data off to mashr
prep = list(Bhat = data.matrix(betas), Shat = data.matrix(beta.ses))
mashrData = mash_set_data(prep$Bhat, prep$Shat)

#1.5)
#account for correlation among samples
sample_corr = estimate_null_correlation_simple(mashrData)
mashrData = mash_update_data(mashrData, sample_corr)

#2)get the covariance matricies between tissue types
#for this specific section, look at this vignette
#https://stephenslab.github.io/mashr/articles/intro_mash_dd.html

#first run ashr to produce a model for the experimental results
#on a per tissue basis.
#From this, pull out any effect with a result significant in any tissue
#This has the effect of increasing the number of tests identified
#as significant
#e.g. for the test data:
#35, 43, 45 sig results in Artery-Aorta, Artery-Tibial, Lung resp. (per original computation)
#(cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Lung/Master.table | grep -P '(chr21\t)' | awk -F"\t" '{if (NR!=1) {print $26}}' | datamash sum 1)
#101, 135, 122 respectively with lfsr < 0.05 after running ashr
#(purrr::map(1:3, function(x) sum(m.1by1$result$lfsr[, x] < 0.05)))
m.1by1 = mash_1by1(mashrData)
strong = get_significant_results(m.1by1,0.05)

#run pca then extreme deconvolution to learn effect patterns from data
#TODO, this should be 5 when not using test data
pca_cov_matrices = cov_pca(mashrData, 2, subset=strong)
ed_cov_matrices = cov_ed(mashrData, pca_cov_matrices, subset=strong) #?? what does this step do?
canonical_cov_matrices = cov_canonical(mashrData)

cov_matrices = c(canonical_co_matrices, ed_cov_matrices)

#3)fit the model
mashrOutput = mash(mashrData, cov_matrices)


