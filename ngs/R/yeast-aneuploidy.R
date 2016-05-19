#!/usr/bin/env Rscript

library(cn.mops)
library(parseArgs)
library(uuutils)
library(ngsutils)
library(rtracklayer)

overview <- function()cat("
Determine whole-chromosome aneuploidy and, if more than one BAM file
is given, determine any 'copy number variations', i.e. regions with
more/fewer reads than expected using cn.mops.

Note: this is yeast specific, with chromosome names like 'chrVI'

Usage:

  yeast-aneuploidy.R  [options]  file1.bam [ file2.bam ... ]

Options:
  --unpaired=TRUE         bam files contain paired-end reads
  --ignore=regions.bed  file with genome regions to ignore

Written by plijnzaad@gmail.com
")

args <- parseArgs(.overview=overview,
                  unpaired=FALSE,
                  ignore=NULL, .allow.rest=TRUE)

bamfiles <- args$.rest
samples <- unname(sapply(bamfiles,
                         function(x)paste(rev(rev(unlist(strsplit(basename(x), "\\.")))[-1]),collapse=".")))

chromos <- 1:16
chromos <- paste0("chr", as.character(as.roman(chromos)))



filter.regions <- filter(counts, iranges, mode=c('any', 'complete')) {
    #counts as returned by cn.mops::getReadCountsFromBAM()
}   #filter.regions

counts <- getReadCountsFromBAM(BAMFiles=bamfiles,
                               mode=ifelse(args$unpaired, "unpaired", "paired"),
                               sampleNames=samples,
                               refSeqName=chromos)

if(!is.null(args$ignore))
  ignore <- import(args$ignore)

counts <- filter.regions(counts, ignore, mode)


chrom.count.stats <- function(bamcounts) {
### bamcounts as returned by cn.mops::getReadCountsFromBAM(a_single_bam_file)
    ## good rule: p < 1e-6 && medianfrac > 0.1
    ## note: this will fail on bamfile containing the yeast RDN1 locus ...
    ## Better filter out anything that is not perfect (bowtie2 mapqual>= 40)
    if( ! any(values(bamcounts)[[1]]>0) )
      stop("only 0's in bamcounts object")
    d <- data.frame(chr=as.factor(seqnames(bamcounts)), counts=as.numeric(values(bamcounts)[[1]]))
    d <- d[ !is.na(d$counts) & d$counts >0 ,]
    grandmean <- mean(d$counts)
    grandmedian <- median(d$counts)
    
    res <- data.frame(pvalue=NA, n=NA, mean=NA, meandiff=NA, median=NA, mediandiff=NA, meanfrac=NA, medianfrac=NA)[0,]
    for(chr in levels(d$chr)) {
        x <- d[d$chr==chr, "counts"]
        n <- length(x)
        t <- wilcox.test(x=(x-grandmedian)) # should use Poisson here? NO
        pval <- p.adjust(t$p.value, method="BH")
        mean <- mean(x)
        meandiff <- mean-grandmean
        meanfrac<- meandiff/grandmean
        md <- median(x)
        mediandiff <- median(x-grandmedian)
        medianfrac <- mediandiff/grandmedian
        res[chr,] <- list(pval, n, mean, meandiff, md, mediandiff, meanfrac,medianfrac)
    }
    res
}                                       # chrom.count.stats

res <- haplocn.mops(counts)
res <- calcIntegerCopyNumbers(res)

#export bed? --bedoutput=
