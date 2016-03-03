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

warning("Running on host ", Sys.getenv("HOSTNAME"), "\n")

### written by plijnzaad@gmail.com
usage <- function(msg) { 
    warning(msg, "\n")
    stop("pattern-content.R \n\
Calculates percentages for fasta file or given genome+location
Usage: env [file=foo.fa | genome=yeast location=chrX:22,000+-100] [ indir=./ ] [ outdir=./] output=filename.{wig,bw,rda} ] pattern=[GC] window=NUMBER  pattern-content.R
")
}

options(verbose=FALSE, stringsAsFactors=FALSE)


indir <- Sys.getenv("indir")
if (indir=="")
  indir <- "."

outdir <- Sys.getenv("outdir")
if (outdir=="")
    outdir <- "."

output <- Sys.getenv("output")
if (output=="")usage("No output specified")

pattern <- Sys.getenv("pattern")
if (pattern=="")usage("Pattern missing")
if (  gregexpr("[^][ACGTUNacgtun]", pattern) > 0)
  usage("Only simple, fixed-length regular expressions without ambiguity codes are allowed")

window <- as.integer(Sys.getenv('window'))
if (is.na(window))usage("no window specified")
  
library(Rsamtools, verbose=FALSE, quietly=TRUE)
library(rtracklayer, verbose=FALSE, quietly=TRUE)
library(ngsutils, verbose=FALSE, quietly=TRUE)
library(zoo, verbose=FALSE, quietly=TRUE)

dna.string <- NULL
chr <- NULL

file <- Sys.getenv('fasta')
file <- paste(indir, file, sep="/")
if ( file.exists(file)) { 
    fa <- open(FaFile(file))
    idx <- scanFaIndex(fa)
    dna <- scanFa(fa,param=idx[1])
    chr <- names(dna)
    dna.string <- as.character(dna) # NOTE: vector of length 1
    stopifnot(length(dna.string) ==1)
}

genome <- Sys.getenv("genome")
if (genome != "") {
    stopifnot(genome=="yeast")
    library(BSgenome.Scerevisiae.UCSC.sacCer3, verbose=FALSE, quietly=TRUE)
    genome <- Scerevisiae
    location <- Sys.getenv('location')
    if(location=="")stop("if genome is specified, location is also required")
    gr <- location2granges(location)
    chr <- as.character(seqnames(gr))
    dna.string <- as.character(genome[[chr]][ranges(gr)])
}
if (is.null(dna.string)) usage("Need either fasta file, or genome + location")

  
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
              seqnames=chr, score=perc)

con <- file(paste0(outdir, "/", output))
if( grepl('.rda$' , con) ) {
    save(file=con, gr)
    warning("Dumping object 'gr'  to ", con, "\nMerge by runnning mergeRda on thesefiles")
} else { 
    warning("Writing output to ", con, "\n")
    export(con=con, object=gr)
    warning("Done. Convert this to BigWig using something like\n\
cat *.wig | wigToBigWig -clip stdin $chromsizes mypattern-71bp.bw\n
or\n
bigWigMerge *.bw out.bedGraph; bedGraphToBigWig out.bedGraph chrom.sizes out.bw\n")
}
sessionInfo()
