#!/usr/bin/env Rscript
##
## Create distribution plot of insert sizes, which are read from the
## tab-delimited file. The sizes can be extracted like
##
##    samtools view -s42.1 foo.bam | awk '!/^@/ && $9 > 0{print $9}' > foo-10perc.insertsizes
## 
## or using e.g. phtools/ngs/sam-insertsizes.pl
##
## issues:
##  - legend quickly unreadable
##
## Note: added to bitbucket.org/princessmaximacenter/rnafusion.git
## on Fri Aug  3 17:56:15 CEST 2018

options(stringsAsFactors = FALSE)
library(optparse)
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(uuutils))

invocation <- paste(commandArgs(), collapse=" ")

usage <- "Usage: %prog [ ] [options] FILENAME(s)

The files should be tab-delimited.  If a filename argument looks like
realfilename.insertlen=COLOR, that color is used; otherwise, rainbow
colors are used."

formalargs <- list(
  make_option("--title", default="insert size distribution", meta='STRING', help='main title of the plot [default: %default]'),
  make_option("--out", default='insert-size-distribution.pdf', meta='FILE', help='name of the output file [default: %default]'),
  make_option("--minlen", default=0L, meta='INTEGER', help='show distribution for inserts no short than this [default: %default]'),
  make_option("--maxlen", default=Inf, meta='INTEGER', help='show distribution for inserts no longer than this [default: %default]'),
  make_option("--add", default=0L, meta='INTEGER', help='adjust insert sizes (e.g., when mapping was done in trimmed way) [default: %default]'),
  make_option("--column", default=1L, meta='INTEGER', help='to find the lengths, use this column in input file [default: %default]'),
  make_option("--scales", default="1", meta='STRING', help='comma-separated list of magnifications for y-axis [default: %default]
                  (This disables log-scaling)'),
  make_option("--no_log", default=FALSE, action='store_true', help='do not use logarithmic y-scale')
  )

epilogue <- ''
parser <- OptionParser(usage=usage,
                     option_list=formalargs,
                       epilogue=epilogue)
p <- parse_args(parser, positional_arguments=TRUE)
opts <- p$options
args <- p$args
rm(p)

opts$log <- !opts$no_log

if( any(grepl(",", opts$scales)) && opts$log) { 
    warning("*** multiscale display is pointless with --log=TRUE, ignored ***")
    opts$no_log <- TRUE
    opts$log <- FALSE
}

## show complete information
print(R.version)
print(Sys.time())
print(invocation)

print("Arguments:\n")
print("\n")

libraries <- c()

## load them, and show where they came from:
for (lib in libraries) { 
  library(lib, character.only=TRUE)
  cat("loaded ", lib, " from ", path.package(lib), "\n")
}

files <- args

if(FALSE) {                             #debugging

    setwd("/Users/philip/tmp/insertsizes/h10000")

    opts <- list()
    opts$title='foo'
    opts$add=0
    opts$minlen=0
    opts$maxlen=1000
    opts$column=1
    opts$out='test.pdf'
    opts$scales="1"                     #must be string
    opts$files <- list.files(pattern="*insert*")
    files <- list.files(pattern="*insert*")
    opts$log <- TRUE

}

log <- opts$log

if (length(files)==0)
  stop("\nNo files.\n\n Run with -h for help")

estimate.mode <- function(x) {          #from uuutils
    d <- density(x)
    d$x[which.max(d$y)]
}                                       #estimate.mode


file.color <- function(file) {
    re <- "^([^=]+)=([^=]+)$"
    if(grepl(re, file, perl=TRUE)) {
        return(unlist(strsplit(file, "=")))
    }
    c(file,"")
}                                       #file.color

all.data <- list()
cum.data <- list()                      #cumulative

colors <- list()
summaries <- list()

for(file in files) {
    fc <- file.color(file)        #pair of color, file
    file <- fc[1]
    color <- fc[2]
    name <- sub("^(.*)\\.[^.]*$", "\\1", file) #strip last extension
    x <- read.table(file)
    if (opts$column > ncol(x) )
      stop("asking for non-existing column")
    data <- x[[opts$column]]
    if (!is.integer(data))
      stop("column should only contain integers")
    if (any(data<0))
      stop("data contains negative values")

    h <- hist(data, nclass=10000, plot=FALSE) # hope this is enough for cumulative
    h$x <- h$mids - 0.5*mean(diff(h$breaks)) #the 'ends' of each bin
    h$y <- cumsum(h$counts)
    h$y <- h$y/max(h$y)           #fraction
    
    cum.data[[name]] <- h

    ## early size selection to keep maximal detail in the density plots
    data <- data[ data >= opts$minlen ]
    data <- data[ data <= opts$maxlen ]
    data <- data + opts$add
    s <- summary(data)
    summaries[[name]] <- c(mean=mean(data), median=median(data), max=max(data), mode=estimate.mode(data))
    all.data[[name]] <- data
    colors[[name]] <- color
}                                       #for files

col <- rainbow(length(all.data), v=0.8)
ncurves <- length(all.data)
lty <- rep(1, ncurves)
if(ncurves>20) { 
    ncolors <- ceiling(ncurves/2)
    col <- rainbow(ncolors, v=0.8)
    col <- as.vector(matrix(rep(col,2), byrow=TRUE, nrow=2))
    lty <- rep(c(1,5), ncolors) # alternate between drawn and dashed, saving on colors
}

names(col) <- names(all.data)
names(lty) <- names(all.data)
for(n in names(colors)) {
    if (colors[[n]] != "") { 
      col[n] <- colors[[n]]
      lty[n] <- 1
  }
}

if (log) {
    densities <- lapply(all.data, function(d){ hist(d, nclass=1000, plot=FALSE)})
    maxy <- max(unlist(lapply(densities,function(d)max(d$counts))))
} else { 
    densities <- lapply(all.data, density)
    maxy <- max(unlist(lapply(densities,function(d)max(d$y))))
}

out <- sub("\\.pdf$","", opts$out)
out <- paste0(out, ".pdf")
title <- opts$title

warning("Creating file ", out)
pdf(file = out, title = title, useDingbats = FALSE, width = 11.7, 
    height = 8.3)

scales <- as.integer(unlist(strsplit(opts$scales, ",")))
stopifnot(length(scales)>0 && sum(is.na(scales))==0)

par(mfrow=c(length(scales), 1), mar=c(4,5,2,4))

minx <- opts$minlen
maxx <- opts$maxlen
if(maxx==Inf)
   maxx <- max(unlist(lapply(all.data,max)))

for(scale in scales) {
    this.maxy <- ifelse(log,log10(maxy), maxy/scale)
    plot(type="n", xlab="",ylab= "",
         x=c(minx, maxx), y=c(0, this.maxy),
         xaxt="n", yaxt="n")
    axis(side=1,at=seq(minx, maxx, 50), labels=TRUE)
    axis(side=1,at=seq(minx, maxx, 10), labels=FALSE, tcl=-0.25) #minor ticks
    abline(v=seq(minx,maxx,100), col="lightgrey")
    mtext("insert size", side=1, line=2.2)
    
    if(log) {
        title(ylab=sprintf("reads"))
        yticks <- 10^(0:6)
        abline(h=1:6, col="lightgrey")
        ylabs <- parse(text=paste("10^", 0:6, sep=""))
        axis(side=2,at=log10(yticks),labels=ylabs, las=1)

    } else {
        title(ylab=sprintf("density x %.0f", scale))
        axis(side=2,labels=FALSE) # labels meaningless since they are densities
    }        

    tcksz <- -0.02
    axis(side=4, tck=tcksz,
         at=seq(0,this.maxy, length.out=11),
         labels=seq(0,100,length.out=11),
         las=1)
    axis(side=4, tck=tcksz/2,              #minor ticks
         at=seq(0,this.maxy, length.out=21),
         labels=NA)
    axis(side=4, tck=tcksz/5,              #minorest ticks
         at=seq(0,this.maxy, length.out=101),
         labels=NA)
    mtext("cumulative percentage", side=4, line=2.2)

    
    if (abs(scale - 1) < 1e-6) {        #put legend in corner of the topmost plot
        title(main=title)
        legend(x="topright", legend=sprintf("%s mean=%.0f med=%.0f mode=%.0f",
                               names(summaries),
                               sapply(summaries,function(x)x['mean']),
                               sapply(summaries,function(x)x['median']),
                               sapply(summaries,function(x)x['mode'])),
               lty=lty[names(all.data)], col=col[names(all.data)])
    }

    for(sample in names(all.data)) {
        if(log) {
            x <- densities[[sample]]$mids
            y <- log10(densities[[sample]]$counts) # real counts (unless truncated)
            nonzero <- densities[[sample]]$counts >0
            lines(lwd=1,
                  x=x[nonzero],
                  y=y[nonzero],
                  col=col[sample],
                  lty=lty[sample])
        } else  {
            d <- densities[[sample]]
            d$y <-  d$y * scale
            lines(d, col=col[sample], lty=lty[sample])
        }
        ## cumulative plots (only  when scale is 1)
        if(abs(scale-1) < 1e-6) { 
            with(cum.data[[sample]],
                 lines(type='s',
                       lwd=1,
                       x=x,
                       y=y*this.maxy, 
                       col=col[sample],
                       lty=lty[sample]))
        }
    }                                   #for samples
}                                       #for scales

dev.off()

warning("Succesfully completed\n")
print("Succesfully completed\n")
sessionInfo()
