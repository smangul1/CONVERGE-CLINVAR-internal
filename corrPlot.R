#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}



#data=read.csv("/u/home/s/serghei/project/CONVERGE/CONVERGE_vs_EXAC_freq/chr1_exacNFE_converge.txt")

data=read.csv(args[1])





png(file = args[2])


plot(data$freqEXAC, data$freqCONVERGE) 
abline(lm(data$freqEXAC~data$freqCONVERGE))
title(args[1])

dev.off()

