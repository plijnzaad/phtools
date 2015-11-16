#!/bin/env Rscript

stop(dput(commandArgs(TRUE)))

## options(verbose=TRUE)

overview <- function()cat(file=stderr(),
                          "Takes all input file (in .bw format), and writes out the sum of their coverages. The intended use is normalization")

library(parseArgs)
args <- parseArgs( .overview=overview,
                  .allow.rest=TRUE,
                  output=NA)

invocation <- paste(commandArgs(), collapse=" ")

## show complete information
print(R.version)
print(Sys.time())
print(invocation)

## libraries to be loaded:
libraries <- c("phutils", "ngsutils", "NCIS", "rtracklayer", "IRanges", "GenomicRanges")

## load them, and show where they came from:
for (lib in libraries) { 
  library(lib, character.only=TRUE)
  cat("loaded ", lib, " from ", path.package(lib), "\n")
}

## options(error=traceback.dump.quit, show.error.messages=T)


print("Arguments:\n")
print(args)
print("\n")

files <- args$.rest
stopifnot(all(file.exists(files)))
stopifnot(length(files)>1)

library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)


## some converting and checking:

stopifnot(all(levels(seqnames(sample1))  == levels(seqnames(sample2))))
stopifnot(all(seqlengths(sample1)==seqlengths(sample2)))

sum <- RleList()

for (cov in covs)  {
  for(chr in names(covs)) { 
    sum[[chr]] <- sum[[chr]]+cov[[chr]]
  }
}
export.bw(sum, con=output)
