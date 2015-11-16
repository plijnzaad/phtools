#!/bin/env Rscript
## 
## Get expression levels given bigwig files (separate for fwd and rev files). Needs work!
## 
## Usage:
##   env fwd=test-fwd.bw rev=test-rev.bw [OR: both=test.bw ] out=test-expr.txt   genome=../cerevisiae-nofasta.gff  bw2expr.R
##

warning("this is not quite ready ... see also bam2count.R")

options(error=traceback.dump.quit, show.error.messages=TRUE)

zero.fill2 <- function(track) { 
  ## make a GRange covering the whole genome, including all the zeros.
  g <- gaps(track)
  values(g)$score <- 0
  sort(c(g, track))
}                                       #zero.fill

get.expression <- function(track, transcripts, FUN=mean, chr=NULL) {
  stopifnot(is(track, "GRanges"))
  stopifnot(is(transcripts, "GRanges"))

  expr <- rep(NA, length(transcripts))
  names(expr) <- names(transcripts)
  counts <- rep(NA, length(transcripts))
  names(counts) <- names(transcripts)
    
  chrs <- if (is.null(chr)) levels(seqnames(track)) else chr
  for (chr in chrs) {
    cat("chr: ", chr, "\n")
    this.strand <- transcripts[seqnames(transcripts)==chr]
    tr <- zero.fill( track[ seqnames(track) == chr  ] )
    scores <- Rle(score(tr), width(ranges(tr)))
    if (length(this.strand)>0) {        #i.e. unless it's the mitochondrial chr
      for(i in 1:length(this.strand) ) {
        ## cat("i: ", i, "\n")
        txpt <- this.strand[i]
        expr[ names(txpt) ] <- FUN(scores[ ranges(txpt) ])
        counts[ names(txpt) ] <- sum(scores[ ranges(txpt) ])
      }
    }
  }
  list(e=expr, c=counts)
}                                       #get.expression

get.expression2 <- function(track, transcripts, FUN=mean) {
  stop("prolly does not work")
## unclear if we can use get.expression2 (i.e. not looping over chromosomes). 
## half of the fwd and half of rev genes are exactly identical, the other half is 
## wrong (in that one of the runs has a zero, the other a non-zero value)
  stopifnot(is(track, "GRanges"))
  stopifnot(is(transcripts, "GRanges"))

  expr <- rep(NA, length(transcripts))
  names(expr) <- names(transcripts)
    
  tr <- zero.fill2(track)
  scores <- Rle(score(tr), width(ranges(tr)))
  if (length(transcripts)>0) {        #i.e. unless it's the mitochondrial chr
      for(i in 1:length(transcripts) ) {
        cat("i: ", i, "\n")
        txpt <- transcripts[i]
        expr[ names(txpt) ] <- FUN(scores[ ranges(txpt) ])
      }
  }
  export.list(expr, counts)
}                                       #get.expression2

genome <- Sys.getenv("genome")
if( ! grepl("\\.gff$", genome))
  stop("no GFF file given for genome\n")
### file <- paste(dir, "cerevisiae-nofasta.gff", sep="/")
genome <- ngsutils::read.sgd.features(file=genome)

expr <- rep(NA, length(genome))
names(expr) <- names(genome)

counts <- rep(NA, length(genome))
names(counts) <- names(genome)

strand.symbols <- c(fwd="+", rev="-", both=NA)

for(strand in c("fwd", "rev", "both")) {
  strand.symbol <- strand.symbols[ strand ]
  track <- Sys.getenv(strand)
  track <- import(track, asRangedData=FALSE)
  use <- if (is.na(strand.symbol)) TRUE else strand(genome) == strand.symbol
  l <- get.expression(track,  transcripts=genome[ use ]) 
  expr [ names(l$e) ] <- l$e
  counts [ names(l$c) ] <- l$c
}

file <- Sys.getenv("out")
if (file=="") {                            #invent name based on that of forward file
  file <- Sys.getenv("fwd")
  if(file=="")
    file <- Sys.getenv("both")
  file <- gsub("(.*)[-:._](fw(d?)|both)\\.bw", "\\1-expr.txt",
               file, ignore.case=TRUE, perl=TRUE)
}

warning("writing to ", file, "\n")
write.tab(file=file, expr[order(names(expr))])


if(TRUE) { 
  file <- "counts.txt"
  write.tab(file=file, counts[order(names(counts))])
  
  file <- "bw2expr.rda"
  warning("dumping to ",file,"\n")
  save.image(file=file)
}

if(FALSE) { 
#### MUCH faster:
  library(Rsubread)

### Rsubread > 1.10 !!!
  annot <- data.frame(GeneID=genome$ID,
                      Chr=as.character(seqnames(genome)),
                      Start=as.integer(start(genome)),
                      End=as.integer(end(genome)),
                      Strand=as.character(strand(genome)))

  counts <- featureCounts(files=file, file.type="BAM",
                          genome=NULL,
                          annot=annot
                          )
}
