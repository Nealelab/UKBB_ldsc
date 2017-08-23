#! /usr/bin/env Rscript

infile <- commandArgs(T)[1]
outname <- commandArgs(T)[2]

dat <- read.table(infile, header=T, stringsAsFactors=F, sep='\t', quote="")
save(dat, file=outname)

summary(dat)

# eof
