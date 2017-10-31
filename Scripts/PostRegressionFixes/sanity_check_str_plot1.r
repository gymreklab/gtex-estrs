#!/usr/bin/env Rscript

## Get mean depth per call from vcf
#system("vcftools --site-mean-depth --gzvcf /storage/szfeupe/Runs/GTEx_estr/HipSTR_OUT/str_calls_chr1.vcf.gz --out Depth")
#system("vcftools --site-depth --gzvcf /storage/szfeupe/Runs/GTEx_estr/HipSTR_OUT/STR_Chr1.filter.vcf.gz --out Depth")
 
## Get distribution of step size for out-of-frame changes and repeat number
#system("awk '{print $4"\t"$10"\t"$5"\t"$6}' /storage/szfeupe/Runs/GTEx_estr/HipSTR_OUT/stutter_models_chr1.txt > step_by_repeat")

## Get the data from depth files into data frames
#call_depth = read.table("Depth.ldepth.mean", sep="\t", header=TRUE)
#call_depth=na.omit(call_depth)
#calls_depth = read.table("Depth.ldepth", sep="\t", header=TRUE)
#calls_depth=na.omit(calls_depth)

## Get step distribution from file and add header
sd = read.table("step_by_repeat", sep="\t", header=FALSE)
names(sd)[1] <- "Step_Dist"
names(sd)[2] <- "Repeats"
names(sd)[3] <- "IDOWN"
names(sd)[4] <- "IUP"

## Aggregate the average distribution by the repeat number
step_dist = aggregate(.~Repeats, data=sd, mean)

## Sanity check of data frames
#head(call_depth)
#head(calls_depth)
head(step_dist)

## Get the plots into standard Rplots.pdf
#plot(density(call_depth$MEAN_DEPTH),  xlim=c(0, 120))
#polygon(density(call_depth$MEAN_DEPTH), col='red', border='blue')
#plot(density(calls_depth$SUM_DEPTH),  xlim=c(0, 1800))

barplot(step_dist$Step_Dist, names.arg=step_dist$Repeats,ylim=c(0, 1.2), xlab="Length of repeats", ylab="Step Distribution")

barplot(step_dist$IDOWN, names.arg=step_dist$Repeats, ylim=c(0,0.15), xlab="Length of repeats", ylab="IDOWN")

barplot(step_dist$IUP, names.arg=step_dist$Repeats, ylim=c(0,0.10), xlab="Length of repeats", ylab="IUP")

#####################Summing IUP and IDOWN
stut_conf = step_dist$IUP + step_dist$IDOWN
barplot(stut_conf, names.arg=step_dist$Repeats, ylim=c(0,0.25), xlab="Size of repeat units", ylab="In frame changes")

plot(step_dist$IUP)

q()

