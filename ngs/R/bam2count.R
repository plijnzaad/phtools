#!/bin/env Rscript
## 
## Get expression levels given bigwig files (separate for fwd and rev files). Needs work!
## 
Usage <- "env [ sam=test.sam | bam=test.bam ] out=test-counts.txt genome=../cerevisiae.gff  bam2count.R"

## l <- commandArgs()
## args <- l[ -(1:grep("^--args$", l)) ]
## if(length(args)==0)
##   stop(Usage)

library(phutils)
library(ngsutils)
library(Rsubread)


stopifnot( packageVersion("Rsubread") > "1.10.4")

file <- Sys.getenv('sam')
if(file =="") {
  file <- Sys.getenv('bam')
  if (file=="")
    stop("Need sam or bam file")
  file.type <- "BAM"
} else {
  file.type <- "SAM"
}

out <- Sys.getenv("out")
stopifnot(out!="")

genome <- Sys.getenv("genome")
if( ! grepl("\\.gff3?$", genome, perl=TRUE))
  stop("no GFF file given for genome\n")

options(error=traceback.dump.quit, show.error.messages=TRUE)


genome <- ngsutils::read.sgd.features(file=genome)

annot <- data.frame(GeneID=genome$ID,
                    Chr=as.character(seqnames(genome)),
                    Start=as.integer(start(genome)),
                    End=as.integer(end(genome)),
                    Strand=as.character(strand(genome)))

counts <- featureCounts(files=file, file.type=file.type, genome=NULL, annot=annot)

write.tab(file=out, data.frame(counts=drop(counts$counts)))
