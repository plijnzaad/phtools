#!/bin/env Rscript

## options(verbose=TRUE)

invocation <- paste(commandArgs(), collapse=" ")

## show complete information
print(R.version)
print(Sys.time())
print(invocation)


library(parseArgs)

overview <- function()cat(file=stderr(),
"normalize.R tries to normalize chip data using NCIS to find the
normalization factor between the ratio between background levels in the
chip and the input data (which should be in .bed format). dir must be
specified; this is the working directory where the input files are to be
found, and (unless outdir is given), where output is written
to. Chromsizes is a file containing chromosome sizes. If e.g.
--chr=chrV,chrVI (comma-separated) is given, only those chromosomes are
read from the input.  --shift=N specifies whether and by how much reads
should be shifted to center the reads (if at all). fraglen is the
initial binwidth. The output is in bigwig format (.bw).  The outputs are
outname-chip-norm.bw for the normalized chip sample (it's scaled to
input, so that one is not output), outname-subtract.bw and
outname-logratio.bw.  If shift was given, outname-chip-cen.bw and
outname-ref-cen.bw are written containing the shifted peaks.")

## parse the command args. NA means mandatory, NULL means optional, rest gets defaults

args <- parseArgs(dir=NA,               # working directory (where input and chip are taken from)
                  outdir=".",          # default: same as dir
                  input=NA,
                  chip=NA,
                  chromsizes=NA,      # tab-delim file with chrName\t\d+
                  outname=NA,
                  chr="all",
                  fraglen=200,
                  shift=0,           #by how much to shift
                  keep_bindata=FALSE,
                  .overview=overview
                  )

print("Arguments:\n")
print(args)
print("\n")

## libraries to be loaded:
libraries <- c("ngsutils", "NCIS", "rtracklayer", "IRanges", "GenomicRanges")

## load them, and show where they came from:
for (lib in libraries) { 
  library(lib, character.only=TRUE)
  cat("loaded ", lib, " from ", path.package(lib), "\n")
}

options(error=traceback.dump.quit, show.error.messages=T)


## some converting and checking:

stopifnot(file_test('-d', args$dir))
setwd(args$dir)

if(is.null(args$outdir))
  args$outdir <- args$dir

stopifnot(file_test('-d', args$outdir))

stopifnot(all(file.exists(c(args$chip, args$input, args$chromsizes))))

if (!is.null(args$shift)) {
  shift <- as.numeric(args$shift)
  stopifnot(is.finite(shift))
  args$shift <- shift
}

{ chromsizes <- args$chromsizes
  stopifnot( nchar(chromsizes)>0 && file.exists(chromsizes) )
  chromsizes <- read.chromSizes(chromsizes)
  stopifnot(all(is.finite(chromsizes)))
  args$chromsizes <- chromsizes
}

if (!is.null(args$chr)) { 
  ## select one or more chromosomes by reading from a pipe running egrep

  args$chr <- unlist(strsplit(args$chr, ","))

  .grep.chr <- function(file, chr) {    
##    tab.re <- '[^-A-Za-z0-9:._/;,]'     # too tricky quoting it
    tab.re <- '[[:space:]]'             # extended character class, hope for the best
    regexp <- sprintf('^(%s)%s[0-9]+' ,
                      paste(chr, sep="", collapse="|"),
                      tab.re)
    pipe <- sprintf(sprintf("egrep '%s' < %s", regexp, file))
##    cat(pipe)
    pipe(pipe, open="r")
  }

  ## pre-select chromosomes here, for speed and more importantly, memory usage
  for (file in c("input", "chip") )
    args[[file]] <- .grep.chr(args[[file]], args$chr)

  print("Preselecting chromosome(s) " ,
          paste(args$chr, collapse=" "),
          " using pipe", "\n")
  args$chr <- NULL
  args$chromsizes
}

args$keep.bindata <- as.logical(args$keep_bindata)


{ fraglen <- as.numeric(args$fraglen)
  stopifnot(is.finite(fraglen))
  args$fraglen <- fraglen 
}

## information about all the arguments:
cat("\n")
print(args)
cat("\n\n")

testing <- as.logical(Sys.getenv("testing"))
if(is.na(testing))testing <- FALSE

chr.len.vec <- args$chromsizes
if (!is.null(args$chr))
  chr.len.vec <- chr.len.vec[args$chr]

if(any(is.na(chr.len.vec)))
  stop("normalize.R: one or more unknown chromosomes among ", paste(args$chr, collapse=" "))

ncis.result <- NCIS(chip.data=args$chip, 
                    input.data=args$input,
                    data.type="BED",
                    frag.len=args$fraglen,
                    chr.vec=args$chr,
                    chr.len.vec=chr.len.vec,
                    keep.bindata=args$keep.bindata)

if(args$keep.bindata)
  save(args, ncis.result, file=sprintf("%s/%s.rda", args$outdir,args$outname))

print("Normalization results:\n\n")
print(ncis.result)
print("\n")

norm.factor <- if (testing) 42 else ncis.result$est

do.export <- function(coverage, dir, name)
  export.bw(coverage, con=sprintf("%s/%s.bw",dir,name))


chip <- ngsutils::get.coverage(file=args$chip, chromsizes=args$chromsizes, shift=args$shift)

if(!is.null(args$shift)) {
  name <- sprintf("%s-chip-cen", args$outname)
  do.export(chip, dir=args$outdir, name)
}

input <- ngsutils::get.coverage(file=args$input, chromsizes=args$chromsizes, shift=args$shift)
if(!is.null(args$shift)) {
  name <- sprintf("%s-ref-cen", args$outname)
  do.export(input, dir=args$outdir, name)
}

.roundit <- function(x)signif(x,digits=5)

norm <- RleList()
for(chr in names(chip))
  norm[[chr]] <- .roundit(chip[[chr]]/norm.factor)
name <- sprintf("%s-chip-norm", args$outname)
do.export(norm, dir=args$outdir, name)


### subtract:
subtract <- RleList()
for(chr in names(chip))
  subtract[[chr]] <- .roundit(chip[[chr]]/norm.factor - input[[chr]])
name <- sprintf("%s-subtract", args$outname)
do.export(subtract, dir=args$outdir, name)

### log ratio:
lograt <- RleList()
for(chr in names(chip))
  lograt[[chr]] <- .roundit(log2(((chip[[chr]]/norm.factor)+1) / (input[[chr]] + 1)))

name <- sprintf("%s-lograt", args$outname)
do.export(lograt, dir=args$outdir, name)

warning("Succesfully completed\n")
print("Succesfully completed\n")
