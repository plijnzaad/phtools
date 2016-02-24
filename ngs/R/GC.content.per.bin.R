#!/usr/bin/env Rscript

### written by plijnzaad@gmail.com

binsize <- as.integer(Sys.getenv("binsize"))
if(is.na(binsize))
  stop("GC.concent.per.bin.R
Calculates GC-content per non-overlapping window along genome, and dumps it as a
GRanges object. Note: genome is assumed to not contain any N's

Usage: env [chr=chrX]  binsize=100   GC.concent.per.bin.R
")

chr.name <- Sys.getenv("chr")

options (stringsAsFactors = FALSE)
library(rtracklayer)
library(ngsutils)
library(uuutils)
## library(GenomeInfoDb)                   #needed for SortSeqlevels

main <- function() {
    library(BSgenome.Scerevisiae.UCSC.sacCer3)
    yeast <- Scerevisiae # defined in this library

    if(chr.name != "")
      rda.file <- sprintf("sacCer3-GCcontent,%s,binsize=%d.rda", chr.name, binsize)
    else
      rda.file <- sprintf("sacCer3-GCcontent,binsize=%d.rda", binsize)
    warning(sprintf("Will dump yeast GC-bins as object 'GC.content' to file %s", rda.file))
    GC.content <- GC.content.of.genome(yeast, binsize=binsize,chr.name=chr.name)
    save(file=rda.file, GC.content)
}

.get.bins <- function(binsize, maxlen) {
    starts <- seq(1, maxlen,by=binsize)
    ends <- starts+binsize-1
    ends[length(ends)] <- maxlen
    cbind(starts,ends)
}

GC.content.per.bin <- function(dna, binsize=50, chr.name, seqinfo) {
    ## bins do not overlap! Assumes there is only A, C, G, and T (no N's)
    ## latter
    stopifnot(is(dna, "DNAString"))
    len <- length(dna)
    bins <- .get.bins(binsize, len)
    r <- apply(bins, 1, function(row){
        letterFrequency(dna[row[1]:row[2]], letters=c("C","G"))})
    r <- apply(r, 2,sum)/binsize
    r [length(r)] <- r [length(r)]*binsize/(len %% binsize) # correct last bin
### r
    GRanges(seqinfo=seqinfo,
            seqnames=chr.name,
            ranges=IRanges(start=bins[,1], end=bins[,2]),
            strand="*",
            score=r)
}                                       #GC.content.per.bin

GC.content.of.genome <- function(genome, binsize=50, chr.names=NULL) { 
    ## bins do not overlap. Assumes there is only A, C, G, and T (no N's)
    stopifnot(is(genome, "BSgenome"))
    result <- list()
    seqinfo <- Seqinfo(seqnames=seqnames(genome), seqlengths=seqlengths(genome), genome="Saccharomyces cerevisiae")
    granges <- GRanges(seqinfo=seqinfo)
    if(is.null(chr.names) || chr.names=="")
      chr.names <- names(genome)
    for(chr in chr.names) {
        gr <- GC.content.per.bin(dna=genome[[chr]],
                                 binsize=binsize,
                                 chr.name=chr,
                                 seqinfo=seqinfo)
        granges <- c(gr,granges)
    }
    ##  do.call("c", sapply(something))  does not work
    ## could have used  BSgenomeViews(genome, gr) here ...
    granges
}                                       #GC.content.of.genome

main()

if(FALSE) { 
### following not used directly, but useful for combining results from previous invocations
### has moved to mergeRData.R
}
