#!/usr/bin/env Rscript

## Create distribution plot of insert sizes, which are read from the
## tab-delimited file(s) (input also be just one column I think).
## 
## The
## files can typically extracted using phtools/ngs/sam-insertlengths.pl
##
## issues:
##  - legend quickly unreadable
##  
##

library(parseArgs)
library(gplots)

invocation <- paste(commandArgs(), collapse=" ")

overview <- function()cat(file=stderr(),
                         "Usage: \
   insertlen-distro.R [ options ] FILENAMES\
\
\
The files should be tab-delimited.  If a filename argument looks like \
realfilename.insertlen=COLOR, \
that color is used; otherwise, rainbow colors are used.\
\
\
Options:\
\
--title=STRING    title of the plots\
--out=FILE        name of output file\
--minlen=INTEGER  show distribution for insert no short than this \
--maxlen=INTEGER  show distribution for insert no longer than this \
--add=INTEGER     adjust insert lengths (e.g., when mapping was done in trimmed way)\
--column=INTEGER  use this column to find the lengths (default: 1)\
--multiscale=STRING    comma-separated list of magnifications (default is to use one)\
                  (This disables log-scaling) 
--log=BOOLEAN     use logarithmic y-scale (default: TRUE)\

")

## --cumulative option has gone since 5-Apr-2018 14:26:16, as cumulative plot is now included by default

args <- parseArgs(out="insert-length-distro.pdf",
                  title='insert length distribution',
                  add=0L,
                  minlen=0L,
                  maxlen=Inf,
                  multiscale="none",
                  log=TRUE,
                  column=1L,
                  .overview=overview,
                  .allow.rest=TRUE
                  )

if (args$multiscale == 'none') { 
  args$scales <- "1"
} else {
    args$scales <- args$multiscale      #old name of the argument
    if(args$log) { 
        warning("*** multiscale display is pointless with --log=TRUE, ignored ***")
        args$log <- FALSE
    }
}

## show complete information
library(uuutils)
print(R.version)
print(Sys.time())
print(invocation)

print("Arguments:\n")
print(args)
print("\n")

libraries <- c()

## load them, and show where they came from:
for (lib in libraries) { 
  library(lib, character.only=TRUE)
  cat("loaded ", lib, " from ", path.package(lib), "\n")
}

files <- args$.rest

if(FALSE) {                             #debugging

    setwd("/Users/philip/tmp/insertsizes/h10000")

    args <- list()
    args$title='foo'
    args$add=0
    args$minlen=0
    args$maxlen=1000
    args$multiscale='none'
    args$column=1
    args$out='test.pdf'
    args$scales="1"                     #must be string
    args$files <- list.files(pattern="*insert*")
    files <- list.files(pattern="*insert*")
    args$log <- TRUE

}

log <- args$log

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

## files <-
##   c("Hsf1t0-1B.insertlen", "Hsf1t0-1D.insertlen=red", "Hsf1t0-2B.insertlen=blue", "Hsf1t0-2D.insertlen=blue")

colors <- list()
summaries <- list()

for(file in files) {
    fc <- file.color(file)        #pair of color, file
    file <- fc[1]
    color <- fc[2]
    name <- sub("^(.*)\\.[^.]*$", "\\1", file) #strip last extension
    x <- read.table(file)
    if (args$column > ncol(x) )
      stop("asking for non-existing column")
    data <- x[[args$column]]
    if (!is.integer(data))
      stop("column should only contain integers")
    if (any(data<0))
      stop("data contains negative values")

    h <- hist(data, nclass=10000, plot=FALSE) # hope this is enough ...
    h$x <- h$mids - 0.5*mean(diff(h$breaks)) #the 'ends' of each bin
    h$y <- cumsum(h$counts)
    h$y <- h$y/max(h$y)           #fraction
    
    cum.data[[name]] <- h

    ## early size selection to keep maximal detail in the density plots
    data <- data[ data >= args$minlen ]
    data <- data[ data <= args$maxlen ]
    data <- data + args$add
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

out <- sub("\\.pdf$","", args$out)
out <- paste0(out, ".pdf")
title <- args$title

warning("Creating file ", out)
pdf(file = out, title = title, useDingbats = FALSE, width = 11.7, 
    height = 8.3)

scales <- as.integer(unlist(strsplit(args$scales, ",")))
stopifnot(length(scales)>0 && sum(is.na(scales))==0)

par(mfrow=c(length(scales), 1), mar=c(2,5,2,3))

minx <- args$minlen
maxx <- args$maxlen
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
    
    if(log) {
        title(ylab=sprintf("reads"))
        yticks <- 10^(0:6)
        abline(h=1:6, col="lightgrey")
        ylabs <- parse(text=paste("10^", 0:6, sep=""))
        axis(side=2,at=log10(yticks),labels=ylabs, las=1)
        title(xlab="fragment length")
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
            lines(lwd=2,
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
