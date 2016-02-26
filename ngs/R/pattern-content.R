#!/bin/env Rscript
##
## script to determine percentage of occurrence of a regular expression
## per base-pair position, averaged over window size, given a fasta
## file. Arguments are all passed as environment variables, invoke as
##
##   env file=fastafile.fa \
##     window=somenumber outdir=somedirectory \
##     pattern='[AG][AG][AG][AG][AG][CT][CT][CT][CT][CT]' pattern-content.R

## (this, BTW, is Trivonof's R5Y5 nucleosome positioning pattern;
## another wellknown 'pattern' would be '[GC]' for simply the GC-content)
## 
## file: fasta file containing ONE chromosome (this is for speed)
## window: length over which to average 
## pattern: regular expression. The pattern will be centered.
##          Only simple, fixed-length REs. No ambiguity codes allowed.
##          Matching is done case-sensitively.
## indir:  where to find fasta files (default: current dir)
## outdir: where to leave results (default: current dir)
##
## This script makes GC-content.R and dinucleotide-content.R obsolete, they
## have therefore been removed as of svn revision 1115
##
## See also pub/tools/sequence/pattern2gff.pl
##


options(verbose=TRUE)

library(Rsamtools)
library(rtracklayer)
library(zoo)

indir <- Sys.getenv("indir")
if (indir=="")
  indir <- "."

outdir <- Sys.getenv("outdir")
if (outdir=="")
    outdir <- "."

pattern <- Sys.getenv("pattern")

if (  gregexpr("[^][ACGTUNacgtun]", pattern) > 0)
  stop("Only simple, fixed-length regular expressions without ambiguity codes are allowed")

window <- as.integer(Sys.getenv('window'))
stopifnot(!is.na(window) && window > 0)

file <- Sys.getenv('file')
file <- paste(indir, file, sep="/")
stopifnot(file.exists(file))

fa <- open(FaFile(file))
idx <- scanFaIndex(fa)
dna <- scanFa(fa,param=idx[1])

chr.name <- names(dna)

dna.string <- as.character(dna) # NOTE: vector of length 1

stopifnot(length(dna.string) ==1)

p <- gsub('[[][^]]+[]]', 'N', pattern)  # e.g. [acg][acg]gg[gc] -> "NNggN"
pat.length <- nchar(p) # lookahead assertion makes match.length attrib all 0!
pattern <- paste0("(?=", pattern, ")")
## Lookahead assertion, needed to also get the overlapping hits
## thanks to http://stackoverflow.com/a/7879329/442359 

hits <- as.vector(gregexpr(pattern, dna.string, perl=TRUE)[[1]]) +
  floor((pat.length-1)/2) ## adjust to middle of the pattern

seq.length <- nchar(dna.string)

occurrence <- rep(0L, seq.length)
occurrence[ hits ] <- 1

perc <- rollmean(x=occurrence, k=window, fill=0)
perc <- as.integer(0.5 + 100*perc)

gr <- GRanges(ranges=IRanges(start=1:seq.length, width=1), strand='*',
              seqnames=chr.name, score=perc)
con <- file(paste0(outdir, "/", chr.name, ".wig"))
warning("Writing output to ", con, "\n")
export(object=gr, con=con, format="wig")
warning("Done. Convert this to BigWig using something like\n\
cat *.wig | wigToBigWig -clip stdin $chromsizes mypattern-71bp.bw\n")
