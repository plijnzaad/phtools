#!/bin/env Rscript
##
## script to determine percentage of occurrence of a regular expression
## per base-pair position, averaged over window size.
## Try e.g. '[GC]' for simply the GC-content, or
## Trivonof's R5Y5 nucleosome positioning pattern: '[AG][AG][AG][AG][AG][CT][CT][CT][CT][CT]'
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
  
library(Rsamtools, verbose=FALSE)
library(rtracklayer, verbose=FALSE)
library(ngsutils, verbose=FALSE)
library(zoo, verbose=FALSE)

dna.string <- NULL
chr <- NULL

seqinfo <- NULL
seqlengths <- NULL

file <- Sys.getenv('fasta')
if(file!="")
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
    seqinfo <- seqinfo(genome)
    seqlengths <- seqlengths(genome)
    location <- Sys.getenv('location')
    if(location=="")stop("if genome is specified, location is also required")
    gr <- location2granges(location, seqinfo=seqinfo, seqlengths=seqlengths)
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
              seqnames=chr, score=perc, seqlengths=seqlengths, seqinfo=seqinfo)

file <- paste0(outdir, "/", output)
if( grepl('.rda$' , file) ) {
    save(file=file, gr)
    warning("Dumping object 'gr'  to ", file, "\nMerge by runnning mergeRda on thesefiles")
} else { 
    warning("Writing output to ", file, "\n")
    export(con=file, object=gr)
    ### NOTE: if exporting bigwig, needs at least the seqlengths!
    warning("Done. Convert this using something like\n\
cat *.wig | wigToBigWig -clip stdin $chromsizes mypattern-71bp.bw\n
or\n
bigWigMerge *.bw out.bedGraph; bedGraphToBigWig out.bedGraph chrom.sizes out.bw\n")
}
sessionInfo()
