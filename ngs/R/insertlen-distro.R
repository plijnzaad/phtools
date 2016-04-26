#!/bin/env Rscript

## issues:
##  - legend quickly unreadable
##  
##

library(parseArgs)

invocation <- paste(commandArgs(), collapse=" ")

overview <- function()cat(file=stderr(),
                         "Usage: \
   insertlen-distro.R [ options ] FILENAMES\
\
The files should be tab-delimited.  If a filename argument looks like realfilename.insertlen=COLOR, \
that color is used; otherwise, rainbow colors are used.\
\
Options:\
\
--title=STRING    title of the plots\
--out=FILE        name of output file\
--maxlen=INTEGER  show distribution up to this length\
--scale=STRING    comma-separated list of magnifications (default: '1,2,5,10,20,50,100')\
--add=INTEGER     adjust insert lengths (e.g., when mapping was done in trimmed way)\
--column=INTEGER  use this column to find the lengths (default: 1)\
--cumulative=BOOLEAN  Print cumulative plots\
--log=BOOLEAN     plot one logarithmic plot, rather than several linear ones. This overrides the --scale option.\

")

args <- parseArgs(out="insert-length-distro.pdf",
                  title='insert length distribution',
                  cumulative=FALSE,
                  add=0L,
                  maxlen=0L,
                  scales="1,2,5,10,20,50,100",
                  log=FALSE,
                  column=1L,
                  .overview=overview,
                  .allow.rest=TRUE
                  )

if(args$log && args$cumulative)
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
    if (!is.integer(x[[args$column]]))
      stop("column should contain integers")
    all.data[[name]] <- x[[args$column]] + args$add
    colors[[name]] <- color
}

col <- rainbow(length(all.data))
names(col) <- names(all.data)
for(n in names(colors)) {
    if (colors[[n]] != "")
      col[n] <- colors[[n]]
}

if(args$cumulative) {
    ecdfs <- lapply(all.data, ecdf)
    maxy <- 1
} else if (args$log) {
    densities <- lapply(all.data, function(d){hist(d, nclass=1000, plot=FALSE}))
    max(y) <- max(unlist(lapply(densities,function(d)max(d$counts))))
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
if(args$log)
  scales <- 1L

par(mfrow=c(length(scales), 1),
    mar=c(1,5,1,1))

maxx <- args$maxlen
if(maxx==0)
   maxx <- max(unlist(lapply(all.data,max)))

for(scale in scales) { 
    plot(type="n", xlab="",ylab= "",
         x=c(0, maxx), y=c(0, ifelse(args$log,maxy/scale, log10(maxy))),
         xaxt="n", yaxt="n")
    axis(side=1,at=seq(0, maxx, 50), labels=TRUE) 
    axis(side=1,at=seq(0, maxx, 10), labels=FALSE, tcl=-0.25) #minor ticks

    if(args$cumulative)
      axis(side=2,labels=TRUE)
    else if(args$log) {
        title(ylab=sprintf("log10(reads)"))
        yticks <- commafy(as.integer(c(1,2,5) %o% 10^(0:5)))
        axis(side=2,at=log10(yticks),label=yticks)
    } else {
        title(ylab=sprintf("density x %.0f", scale))
        axis(side=2,labels=FALSE)
    }        
    
    if (abs(scale - 1) < 1e-6) {        #put legend in corner of the topmost plot
        title(main=title)
        legend(x="topright", legend=names(all.data),
               lty=1, col=col[names(all.data)])
    }

    for(sample in names(all.data)) {
        if(args$cumulative) {
            f <- ecdfs[[sample]]
            x <- seq(0, maxx, length.out=500)
            y <- f(x)*scale
            lines(x,y, col=col[sample])
        } else if(arg$log) {
            lines(x=densities[[sample]]$mids, y=log10(densities[[sample]]$counts), col=col[sample])
        } else  {
             d <- densities[[sample]]
             d$y <-  d$y * scale
             lines(d, col=col[sample])
         }
    }
}                                       #for scale

dev.off()

warning("Succesfully completed\n")
print("Succesfully completed\n")
sessionInfo()
