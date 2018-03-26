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
--cumulative=BOOLEAN  Print cumulative plots\
--multiscale=STRING    comma-separated list of magnifications (default is to use one logarithmic scale)\

")

args <- parseArgs(out="insert-length-distro.pdf",
                  title='insert length distribution',
                  cumulative=FALSE,
                  add=0L,
                  minlen=0L,
                  maxlen=Inf,
                  multiscale="none",
                  column=1L,
                  .overview=overview,
                  .allow.rest=TRUE
                  )

log <- TRUE
if (args$multiscale == 'none') { 
    args$scales <- "1"
} else {
    log <- FALSE
    args$scales <- args$multiscale      #old name of the argument
}

if(log && args$cumulative)
  stop("logarithmic and cumulative prolly don't work, first check this")

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

file.color <- function(file) {
    re <- "^([^=]+)=([^=]+)$"
    if(grepl(re, file, perl=TRUE)) {
        return(unlist(strsplit(file, "=")))
    }
    c(file,"")
}                                       #file.color

all.data <- list()

## files <-
##   c("Hsf1t0-1B.insertlen", "Hsf1t0-1D.insertlen=red", "Hsf1t0-2B.insertlen=blue", "Hsf1t0-2D.insertlen=blue")

colors <- list()

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
      stop("column should contain integers")
    if (any(data<0))
      stop("data contains negative values")
     # early size selection to keep maximal detail in the density plots
    data <- data[ data >= args$minlen ]
    data <- data[ data <= args$maxlen ]
    all.data[[name]] <- data + args$add
    colors[[name]] <- color
}

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

if(args$cumulative) {
    ecdfs <- lapply(all.data, ecdf)
    maxy <- 1
 } else if (log) {
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
if(log) { 
  pdf(file = out, title = title, useDingbats = FALSE, width = 11.7, 
      height = 8.3)
} else  { 
  pdf(file = out, title = title, useDingbats = FALSE, height = 11.7, 
      width = 8.3)
}

scales <- as.integer(unlist(strsplit(args$scales, ",")))
stopifnot(length(scales)>0 && sum(is.na(scales))==0)

if(log) { 
    scales <- 1L
} else { 
    par(mfrow=c(length(scales), 1), mar=c(1,5,1,1))
}

minx <- args$minlen
maxx <- args$maxlen
if(maxx==Inf)
   maxx <- max(unlist(lapply(all.data,max)))

for(scale in scales) { 
    plot(type="n", xlab="",ylab= "",
         x=c(minx, maxx), y=c(0, ifelse(log,log10(maxy), maxy/scale)),
         xaxt="n", yaxt="n")
    axis(side=1,at=seq(minx, maxx, 50), labels=TRUE)
    axis(side=1,at=seq(minx, maxx, 10), labels=FALSE, tcl=-0.25) #minor ticks
    abline(v=seq(minx,maxx,100), col="lightgrey")
    
    if(args$cumulative) {
        axis(side=2,labels=TRUE)
        abline(v=seq(h,1,0.1), col="lightgrey")
    } else if(log) {
        title(ylab=sprintf("reads"))
        yticks <- 10^(0:6)
        abline(h=1:6, col="lightgrey")
        ylabs <- parse(text=paste("10^", 0:6, sep=""))
        axis(side=2,at=log10(yticks),labels=ylabs, las=1)
        title(xlab="fragment length")
    } else {
        title(ylab=sprintf("density x %.0f", scale))
        axis(side=2,labels=FALSE)
    }        
    
    if (abs(scale - 1) < 1e-6) {        #put legend in corner of the topmost plot
        title(main=title)
        legend(x="topright", legend=names(all.data),
               lty=lty[names(all.data)], col=col[names(all.data)])
    }

    for(sample in names(all.data)) {
        if(args$cumulative) {
            f <- ecdfs[[sample]]
            x <- seq(minx, maxx, length.out=500)
            y <- f(x)*scale
            lines(x,y, col=col[sample], lty=lty[sample])
        } else if(log) {
            lines(x=densities[[sample]]$mids, y=log10(densities[[sample]]$counts), col=col[sample], lty=lty[sample])
        } else  {
             d <- densities[[sample]]
             d$y <-  d$y * scale
             lines(d, col=col[sample], lty=lty[sample])
         }
    }
}                                       #for scale

dev.off()

warning("Succesfully completed\n")
print("Succesfully completed\n")
sessionInfo()
