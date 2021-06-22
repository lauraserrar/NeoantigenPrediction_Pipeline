#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two argument must be supplied (input file, sample).n", call.=FALSE)
}

library(sequenza)

Sys.time()
test <- sequenza.extract(args[1])
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = args[2], out.dir = args[2])

Sys.time()