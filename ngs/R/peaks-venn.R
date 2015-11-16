#!/bin/env Rscript

options(verbose=FALSE)
library(parseArgs, quietly=TRUE, verbose=FALSE)

overview <- function()
    cat(file=stderr(),"

  Usage: peaks-venn.R [ options ] foo.bed bar.bed [ ... ]

  Creates Venn diagram of overlaps amongst peaks in (up to 4) bed files
  The directory and the '_peaks.bed' parts of the file names are ignored for
  brevity. Originally written for MACS2 output, but should also work for
  other .bed files. In addition, the peaks belonging to each of the different
  overlaps amongst the input sets can be written to separate .bed files.

  Options:

  --out=FILE.PDF  Specify output file (default: peaks-venn.pdf)

  --strip_pdir=FALSE  When naming the sets, ignore the directory component
                      of their input files (default: TRUE)

  --peakset_prefix=STRING Write the various overlap sets to STRING-subset.bed
                          (default: '', i.e. don't dump anything).

  --verbose=BOOLEAN    Whether to be verbose (default: FALSE)
")

args <- parseArgs(out="peaks-venn.pdf",
                  strip_dir=TRUE,
                  peakset_prefix="",
                  verbose=FALSE,
                  .allow.rest=TRUE,
                  .overview=overview)

output <- args$out
strip.dir <- args$strip_dir
peakset.prefix <- args$peakset_prefix
files <- args$.rest
options(verbose=args$verbose)

if(length(files) < 2 || length(files) > 4)
    stop("Need more than 1 and fewer than 4 .bed-files")

library(DiffBind, quietly=TRUE, verbose=args$verbose)
stopifnot(packageVersion("DiffBind") == "1.14.99") # 1.14.5 etc. have empty subset bug
library(rtracklayer, verbose=args$verbose)

ids <- gsub("_peaks", "", files)
ids <- gsub(".bed", "", ids)
if(strip.dir)
  ids <- gsub(".*/", "", ids)
cat(ids, "\n")
peaksets <- data.frame(SampleID=ids,
                       Peaks=files,
                       PeakCaller="bed")
dba <- dba(sampleSheet=peaksets, bCorPlot=FALSE)

warning("creating file ", output ,"\n")
pdf(file = output, title = file, useDingbats = FALSE,
    width = 8.3, height = 11.7)

results <- dba.plotVenn(dba, mask=rep(TRUE, nrow(peaksets)),
                        bReturnPeaksets= TRUE )

dev.off()

### reformat peaksets so we can dump it to stdout
peaksets$letter <- LETTERS[1:nrow(peaksets)]
peaksets$file <- peaksets$Peaks
peaksets$size <- unlist(lapply(dba$peaks, nrow))
peaksets$PeakCaller <- NULL
peaksets$Peaks <- NULL
cat("========================================================================\n")
write.table(x=peaksets, file=stdout(), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
cat("========================================================================\n")
for(subset in names(results)) {
    size <- if(is.null(results[[subset]])) 0 else length(results[[subset]])
    cat(subset, "\t", size, "\n")
}
cat("========================================================================\n")

if (peakset.prefix != "")
  for(subset in names(results)) {
      file <- paste0(peakset.prefix, "-", subset, ".bed")
      if(is.null(results[[subset]])) {
          warning("subset ", subset," is empty, created empty file\n")
          file.create(file)
      } else {
          warning("dumping subset to ", file, "\n")
          export(results[[subset]], con=file,format="bed")
      }
  }

if(args$verbose)
  sessionInfo()
