#!/bin/env Rscript-3.2.1

options(verbose=FALSE)
stopifnot( getRversion() >= "3.2.1")
library(parseArgs, quietly=TRUE, verbose=FALSE)

### options and their defaults:
fragmentLen <- 147
width <- 3
percentile <- 30
trim <- 73
region <- ""
pcKeepComp <- 'auto'                    #between 0.01 and 0.03: hardly any influence
verbose <- FALSE

overview <- function()
    cat(file=stderr(),"
Usage: call-nucleosomes.R  --bam=input.bam --out=output.bed [ options ]

   To use standard output, use --out=stdout. You cannot use standard input.

Options:

  --region=VALUE       Genome region (in samtools format) to limit the analysis
                       to. Use whitespace to separate several regions
  --fragmentLen=VALUE  Use this as the nucleosome size in bp. 0 means: estimate
                       from data (very slow). Note that if reads where trimmed at
                       their 5'-end during mapping, fragmentLen also should be less
  --trim=VALUE         Shorten reads to this length, to sharpen the peaks
  --pcKeepComp=VALUE   Use this fraction of frequency components. 'auto' makes
                       nucleR choose this. Rarely needs changing
  --percentile=VALUE   Only peaks with score in greater than this percentile are called
  --width=VALUE        Output peaks with this width
  --verbose=BOOLEAN    Whether to be verbose



")

args <- parseArgs(bam=NA,
                  out=NA,
                  region=region,
                  fragmentLen=fragmentLen,
                  trim=trim,
                  pcKeepComp=pcKeepComp,
                  percentile=percentile,
                  width=width,
                  verbose=verbose,
                  .allow.rest=FALSE,
                  .overview=overview)

library(uuutils,quietly=!verbose)
import.list(args)                       #re-assigns things under same name

options(verbose=verbose)

library(GenomicAlignments, quietly=!verbose)
library(ShortRead, quietly=!verbose)
library(ngsutils, quietly=!verbose)
library(rtracklayer, quietly=!verbose)
library(nucleR, quietly=!verbose)
stopifnot(packageVersion("nucleR") >= "2.0.99")

## options(error=traceback.dump.quit(), show.error.messages=TRUE)

## if(bam == "stdin") ### does not work ...
##   bam <- stdin()

region <- unlist(strsplit(region, " "))

allreads <- read.bam.region(file=bam, region=region)
allreads <- as(allreads, "RangedData")

allpeaks <- GRanges()

debug <- FALSE

for(chr in unique(seqnames(allreads))) {
    ## select, trim and shift
    reads <- processReads(allreads[chr], type="single",
                          fragmentLen=fragmentLen, trim=trim)
    coverage <- coverage.rpm(reads)
    coverage <- filterFFT(coverage[[1]], pcKeepComp=pcKeepComp, showPowerSpec=debug)
    peaks <- peakDetection(coverage, threshold=paste(percentile, "%"),
                           score=TRUE, width=width)
    seqlevels(peaks) <- seqlevels(allreads)
    names(peaks) <- chr                 # (hack)
    gr <- as(peaks, "GRanges")
    allpeaks <- c(allpeaks, gr)
}

if(out=='stdout')
  out <- stdout()

rtracklayer::export(peaks, con=out, format="bed")


