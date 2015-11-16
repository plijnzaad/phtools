#!/bin/env Rscript

## script to determine k-mer frequency per base-pair position of a fasta file.
## Expects k in env.var. 'k' and name of fasta file in env.var. 'chr' (should
## contain one chromosome)

options(verbose=TRUE)

library(Rsamtools)
library(combinat)

seqdir <- Sys.getenv("seqdir")    # where to find sequence files
if(seqdir=="")
  seqdir="."

outdir <- Sys.getenv("outdir")    # where to write output (can be empty)

alphabet <- c("A","C", "G", "T")        #could extend this with e.g. 'N'

all.kmers <- function(alphabet, k) { 
  n.alphabet <- length(alphabet)
  h <- hcube(rep(n.alphabet, k))        # bit like expand.grid
  dim <- dim(h)
  h <- alphabet[ h ]
  dim(h) <- dim
  decoding <- apply(h, 1, paste, collapse="")
  coding <- 1L:(n.alphabet^k)
  names(coding) <- decoding
  coding
}

k <- as.integer(Sys.getenv('k'))
stopifnot(k>1)

coding <- all.kmers(alphabet=alphabet, k=k)
decoding <- names(coding)

chr <- Sys.getenv('chr')
file <- paste(seqdir, chr, sep="/")
stopifnot(file.exists(file))

fa <- open(FaFile(file))  
idx <- scanFaIndex(fa)
dna <- scanFa(fa,param=idx[1])

chr.name <- names(dna)

dna.string <- as.character(dna)
stopifnot(length(dna.string) ==1)
dna.chars <- unlist(strsplit(dna.string, "")) # one letter per array elt

n <- length(dna.chars) - k + 1 

m <- sapply(1:k, function(i)dna.chars[(1:n)+i-1])
shards <-  apply(m, 1, paste, collapse="")

seq <- unname(coding[ shards ] )        # seq of integers, indicating kmer (left-aligned!)

kmers <- list()

kmers$k <- k                           #window size
kmers$chr.name <- chr.name
kmers$kmers <- seq                      #sequence of integers
kmers$coding <- coding                  # coding["GAA"] => 3
kmers$decoding <- decoding              # e.g. decoding[3] => GAA

file <- paste0(outdir,chr, '.rda')
cat("Writing results to ", file, "\n")
save(file=file, kmers)
