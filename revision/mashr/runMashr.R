#!/usr/bin/env Rscript

testrun=FALSE
testrun2Chroms=FALSE
if (testrun) {
	workdir = 'testrun'
} else if (testrun2Chroms) {
	workdir = 'testrun2Chroms'
} else {
	workdir = 'fullrun'
}
indir = paste(workdir, '/input', sep='')
outdir = paste(workdir, '/output', sep='')
intermediate = paste(workdir, '/intermediate', sep='')

loadData = function(indir, intermediate) {
	print('----Load the data into two dataframes----')
	library('purrr')
	#get files to load
	files = list.files(path=indir, pattern="*.table", full.names=TRUE)
	shortFileNames = purrr::map(files, function(name) basename(tools::file_path_sans_ext(name)))

	#load files individually
	#lung.table <- read.table("Lung.table", header=TRUE, colClasses=c('character', 'factor', 'integer', 'numeric', 'numeric'))
	dataFrameList = purrr::map(files, read.table, header=TRUE, colClasses=c('character', 'factor', 'integer', 'numeric', 'numeric'))

	#some rows are duplicated for some reason, so remove them
	dataFrameList = purrr::map(dataFrameList, function(x) x[!duplicated(x), ])

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
	#In the full run, this will be sized: 226562 obs. of  37 variables:
	allData = purrr::reduce(dataFrameList, function(df1, df2) base::merge(df1, df2, by=c('gene', 'chrom', 'str.start')))

	#move gene, str idetnifier info into rownames
	rownames(allData) = purrr::pmap(allData[c(1,2,3)], paste, sep='_')
	allData = allData[-(1:3)]

	#separate betas and beta.ses
	betas = allData[!grepl('beta.se', colnames(allData))]
	beta.ses = allData[grepl('beta.se', colnames(allData))]

	colnames(betas) = purrr::map(colnames(betas), substring, 6)
	colnames(beta.ses) = purrr::map(colnames(beta.ses), substring, 9)

	saveRDS(betas, paste(intermediate, '/betas.rds', sep=''))
	saveRDS(beta.ses, paste(intermediate, '/beta.ses.rds', sep=''))
	return(list(betas, beta.ses))
}

prepMashr = function(betas, beta.ses, intermediate) {
	print('----prep mashr----')
	library(mashr)
	#For overall coding pattern, look at this vignette:
	#https://stephenslab.github.io/mashr/articles/intro_mash.html

	#1)hand the data off to mashr
	prep = list(Bhat = data.matrix(betas), Shat = data.matrix(beta.ses))
	mashrData = mash_set_data(prep$Bhat, prep$Shat)

	#1.5)
	#account for correlation among samples
	sample_corr = estimate_null_correlation_simple(mashrData)
	mashrData = mash_update_data(mashrData, V=sample_corr)
	saveRDS(mashrData, paste(intermediate, '/mashrData.rds', sep=''))
	saveRDS(sample_corr, paste(intermediate, '/sample_corr.rds', sep=''))

	#1.6)
	#should we be doing this: https://stephenslab.github.io/mashr/articles/intro_mashnobaseline.html ?

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
	if (testrun || testrun2Chroms) {
		num_components = 2 #not enough tissues to use 5 components in the test run
	} else {
		num_components = 5
	}
	pca_cov_matrices = cov_pca(mashrData, num_components, subset=strong)
	ed_cov_matrices = cov_ed(mashrData, pca_cov_matrices, subset=strong) #?? what does this step do?
	canonical_cov_matrices = cov_canonical(mashrData)

	cov_matrices = c(canonical_cov_matrices, ed_cov_matrices)
	saveRDS(cov_matrices, paste(intermediate, '/cov_matrices.rds', sep=''))

	#3)fit the model and save its output

	mashrModel = mash(mashrData, cov_matrices, outputlevel = 1)
	saveRDS(mashrModel, paste(intermediate, '/mashrModel.rds', sep=''))	
	fittedG = get_fitted_g(mashrModel)
	saveRDS(fittedG, paste(intermediate, '/fittedG.rds', sep=''))	

	return(list(mashrData, fittedG, sample_corr))
}

runMashr = function(mashrData, fittedG, outdir) {
	print(paste('----run mashr with outdir ', outdir, '----', sep=''))
	library(mashr)

	mashrOutput = mash(mashrData, g=fittedG, fixg=TRUE)
	saveRDS(mashrOutput, paste(outdir, '/mashrOutput.rds', sep=''))

	lfsr = get_lfsr(mashrOutput)
	write.table(lfsr, paste(outdir, '/posterior_lfsr.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	posterior_betas = get_pm(mashrOutput)
	write.table(posterior_betas, paste(outdir, '/posterior_betas.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	posterior_beta_ses = get_psd(mashrOutput)
	write.table(posterior_beta_ses, paste(outdir, '/posterior_beta_ses.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	log10bf = get_log10bf(mashrOutput)
	rownames(log10bf) = rownames(posterior_beta_ses)
	write.table(log10bf, paste(outdir, '/posterior_log10bf.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)

	#4)
	#Create correlation plots heatmaps
	#From this vignette https://stephenslab.github.io/mashr/articles/mash_sampling.html
	#would like to use the sampling based method - should be more conservative, but it uses R instead of C backing and this seems to slow (could try waiting)
	library(corrplot)
	sharing = get_pairwise_sharing(mashrOutput, factor = 0.5)
	png(paste(outdir, '/significantEffectSharing.png', sep=''))
	corrplot(sharing, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
	dev.off()

	sig = get_significant_results(mashrOutput)
	secorr = cor(get_pm(mashrOutput)[sig, ])
	png(paste(outdir, '/significantEffectCorr.png', sep=''))
	corrplot(secorr, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
	dev.off()
}

runMashrChromByChrom = function(betas, beta.ses, sample_corr, fittedG, outdir) {
	print('---running mashr chrom by chrom---')
	chroms = paste("chr", c(1:22, 'X'), sep='')
	for (chrom in chroms) {
		rows = grepl(paste(chrom, '_', sep=''), rownames(betas), fixed=TRUE)
		if (!any(rows)) {
			#in some tests, not all chromosomes will be present
			next
		}
		chromBetas = betas[rows, ]
		chromBeta.ses = beta.ses[rows, ]
		prep = list(Bhat = data.matrix(chromBetas), Shat = data.matrix(chromBeta.ses))
		chromMashrData = mash_set_data(prep$Bhat, prep$Shat, V=sample_corr)
		chromOutdir = paste(outdir, chrom, sep='/')
		dir.create(chromOutdir, showWarnings = FALSE)
		runMashr(chromMashrData, fittedG, chromOutdir)
	}
}

#Step 1: load data
l = loadData(indir, intermediate)
betas = l[[1]]
beta.ses = l[[2]]
#betas = readRDS(paste(intermediate, '/betas.rds', sep=''))
#beta.ses = readRDS(paste(intermediate, '/beta.ses.rds', sep=''))

#Step 2: prep mashr
l = prepMashr(betas, beta.ses, intermediate)
mashrData = l[[1]]
fittedG = l[[2]]
sample_corr = l[[3]]
#mashrData = readRDS(paste(intermediate, '/mashrData.rds', sep=''))
#fittedG = readRDS(paste(intermediate, '/fittedG.rds', sep=''))
#sample_corr = readRDS(paste(intermediate, '/sample_corr.rds', sep=''))

#step 3: run mashr
#runMashr(mashrData, fittedG, outdir)
runMashrChromByChrom(betas, beta.ses, sample_corr, fittedG, outdir)

#TODO : run without the various optimizations to see what changes
#Are we using EZ model or EE model? Mentioned in paper
#The paper initializes extreme deconvolution in a more sophisticated manner than the above
	#Should we do that?


