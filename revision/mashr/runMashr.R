#!/bin/R

#------------------Load the data into one big dataframe----------------
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

allData = purrr::reduce(dataFrameList, function(df1, df2) base::merge(df1, df2, by=c('gene', 'chrom', 'str.start')))

#-------------------------Hand off to mashr------------------------
library(mashr)

#move gene, str idetnifier info into rownames
rownames(allData) = purrr::pmap(allData[c(1,2,3)], paste, sep='_')
allData = allData[-(1:3)]
betas = allData[!grepl('se', colnames(allData))]
beta.ses = allData[grepl('se', colnames(allData))]

#hand the data off to mashr
prep = list(Bhat = data.matrix(betas), Shat = data.matrix(beta.ses))
mashrData = mash_set_data(prep$Bhat, prep$Shat)

#get the covariance matricies
U.c = cov_canonical(mashrData)

#fit the model
m.c = mash(mashrData, U.c)


