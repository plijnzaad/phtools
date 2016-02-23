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
    for(chr in names(genome))
      result[[chr]] <- GC.content.per.bin(genome[[chr]], binsize=binsize)
    result
}                                       #GC.content.of.genome

.get.bins <- function(binsize, maxlen) {
    starts <- seq(1, maxlen,by=binsize)
    ends <- min(starts+binsize-1, maxlen)
    cbind(starts,ends)
}

GC.content.per.bin <- function(dna, binsize=50) {
    ## bins do not overlap! Assumes there is only A, C, G, and T (no N's)
    ## latter
    stopifnot(is(dna, "DNAString"))
    len <- length(dna)
    ### bins <- seq(1, len,by=binsize)
    bins <- .get.bins(binsize, len)
    r <- apply(bins, 1, function(row){
        letterFrequency(dna[row[1]:row[2]], letters=c("C","G"))
    })
    ### r <- sapply(bins, function(bin){
    ###         letterFrequency(dna[bin:(min(bin+binsize-1, len))], letters=c("C","G"))
    ###     })
    r <- apply(r, 2,sum)/binsize
    # correct last bin:
    r [length(r)] <- r [length(r)]*binsize/(len %% binsize)
    r
}                                       #GC.content.per.bin
