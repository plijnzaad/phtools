### Create tracks of GC content per (non-overlapping) window along genome
### written by plijnzaad@gmail.com

options (stringsAsFactors = FALSE)
library(rtracklayer)
library(ngsutils)
library(uuutils)
## library(GenomeInfoDb)                   #needed for SortSeqlevels

library(BSgenome.Scerevisiae.UCSC.sacCer3)
yeast <- Scerevisiae # defined in this library

GC.content.of.genome <- function(genome, binsize=50) { 
    ## bins do not overlap. Assumes there is only A, C, G, and T (no N's)
    stopifnot(is(genome, "BSGenome"))
    result <- list()
    seqinfo <- Seqinfo(seqnames=seqnames(yeast), seqlengths=seqlengths(yeast), genome="Saccharomyces cerevisiae")
    granges <- GRanges(seqinfo=seqinfo)
    for(chr in names(genome)) { 
        gr <- GC.content.per.bin(dna=genome[[chr]],
                                 binsize=binsize,
                                 chr.name=chr,
                                 seqinfo=seqinfo)
        result[[chr]] <- gr
    }
    do.call("c", result)
}                                       #GC.content.of.genome

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
