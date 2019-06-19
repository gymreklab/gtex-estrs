How to run:

# TODO SNPs

./gatherData.sh
nohup R CMD BATCH '--args runval="strs"' runMashr.R /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/mashr.strs.log &
# To follow: tail -f /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/mashr.strs.log



Examining output and other files:
outputDirectory/output - contains tsvs of the four posterior calculations
	(betas, beta_ses, lfsr and log10 of the bayes factors),
	plus individual directories per chromosome
outputDirectory/input - a cleaned copy made of the data the script is using
outputDirectory/intermediate - intermediate R objects generated during
	running runMashr.R - useful for resuming from an intermediate state
	or inspecting what's going on

If you wish to change the output directory:
Change the workdir variable on line 71 of gatherData.sh
Change the workdir variable on line 12 of runMashr.r

If you wish to change the input directory:
Change the tissuedir variable on line 96 of gatherData.sh
Possibly change line 98 of gatherData.sh - either
	pointing the cat command to the correct files,
	or making sure that awk is pulling out the correct
	columns (gene, chrom, str.id, beta, beta_se)

If you wish to work on the extreme deconvolution bug:
Beginning at line 44, here's what I think the code should be
#Notice all=TRUE is uncommented
allData = purrr::reduce(dataFrameList, function(df1, df2) base::merge(df1, df2, by=c('gene', 'chrom', 'str.id'), all=TRUE))     #comment out this line to keep NAs instead of zeros for STRs with results in only some tissues
#by default, missing data is NA, set it to something small and irrelevant
allData[is.na(allData)] = 1e-7

