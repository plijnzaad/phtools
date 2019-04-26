#!/bin/env Rscript
## script to determine k-mer frequency per base-pair position of a fasta
## file.  Expects k in env.var. 'k' and name of fasta file in
## env.var. 'chr' (should contain one chromosome), dumps list of kmers
## found and their (relative) abundances to Rdata file

options(verbose=TRUE)

library(Rsamtools)
library(combinat)

seqdir <- Sys.getenv("seqdir")    # where to find sequence files
if(seqdir=="")
  seqdir="."

outdir <- Sys.getenv("outdir")    # where to write output (can be empty)

alphabet <- c("A", "C", "G", "T")      # could extend this with e.g. 'N'

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

# in matrix form (each row  == one kmer)
kmer.mat <- unlist(strsplit(decoding,""))
kmer.mat <- matrix(kmer.mat, ncol=3, byrow=TRUE)

chr <- Sys.getenv('chr')
file <- paste(seqdir, chr, sep="/")
stopifnot(file.exists(file))

fa <- open(FaFile(file))
idx <- scanFaIndex(fa)
dna <- scanFa(fa,param=idx[1])

chr.name <- names(dna)

dna.string <- as.character(dna)
stopifnot(length(dna.string) ==1)
dna.chars <- unlist(strsplit(dna.string, "")) # vector with one letter per nt

n <- length(dna.chars) - k + 1 

m <- sapply(1:k, function(i)dna.chars[(1:n)+i-1]) # row = length k subsequence
shards <-  apply(m, 1, paste, collapse="") # same but has strings of lenght k

seq <- unname(coding[ shards ] )        # seq of integers, indicating kmer (left-aligned!)


## nucleotide composition:
composition <- table(m[,1])
frac.composition <- composition/sum(composition)

## given this nt composition, what is the likelihood of each kmer 
## occurring by chance:

bg.likelihoods <- apply(kmer.mat, 1, function(row)prod(frac.composition[row]))
names(bg.likelihoods) <- decoding
## NOTE: sum of these is only 1 if nt composition is uniform (which it isn't)
## double check if probabilities sum to 1:
## frac <- rep(0.25 , 4); names(frac) <- names(frac.composition)
## sum(apply(kmer.mat, 1, function(row)prod(frac[row]))) => 1

## scale to 1:
names(bg.likelihoods) <- decoding
bg.likelihoods <- bg.likelihoods/sum(bg.likelihoods)

## occurence per kmer, sorted: 
tab <- table(seq)
names(tab) <- decoding[as.integer(names(tab))]
tab <- sort(tab,decr=TRUE)
tab.frac <- tab/(sum(tab))
tab.perc <- 100*tab.frac

## occurrence relative to random background:
rel.abundance <- sort(tab.frac/bg.likelihoods[names(tab.frac)], decr=TRUE)

kmers <- list()                         # object to contain all the info

kmers$k <- k                            # window size
kmers$chr.name <- chr.name
kmers$coding <- coding                  # coding["GAA"] => 3
kmers$decoding <- decoding              # e.g. decoding[3] => GAA

kmers$frac.composition <- frac.composition #of query sequence
kmers$bg.likelihoods <- bg.likelihoods     # per k-mer

kmers$kmers <- seq                      # sequence of integers pointing to kme
kmers$tab <- tab                        # kmers sorted by occurrence
kmers$tab.frac <- tab.frac              # same but fractional
kmers$rel.abundance <- rel.abundance    # over/under representation per k-mer

# save results as an R data object
file <- paste0(outdir,chr, '.rda')
cat("Writing results to ", file, "\n")
save(file=file, kmers)
