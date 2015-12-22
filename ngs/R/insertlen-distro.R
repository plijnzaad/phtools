#!/bin/env Rscript

## issues:
##  - legend is bit unreadable
##  - tick lengths at the labels must be longer
##  - allow lwd, lty too?
##

library(parseArgs)

invocation <- paste(commandArgs(), collapse=" ")

overview <- function()cat(file=stderr(),
                         "Usage: \
   insertlen-distro.R [ options ] FILENAMES\
\
If a file name argument looks like realfilename.insertlen=COLOR, that color is used; otherwise,\
rainbow colors are used.\
\
Options:\
\
--title     title of the plots\
--out       name of output file\
--maxlen    show distribution up to this length\
--scales    comma-separated list of magnifications\
--add       adjust insert lengths (e.g., when mapping was done in trimmed way)\

")

args <- parseArgs(out="insert-length-distro.pdf",
                  title='insert length distribution',
                  add=0L,
                  maxlen=0L,
                  scales="1,2,5,10,20,50,100",
                  .overview=overview,
                  .allow.rest=TRUE
                  )

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
    stopifnot(ncol(x)==1)
    stopifnot(is.integer(x[[1]]))
    all.data[[name]] <- x[[1]] + args$add
    colors[[name]] <- color
}

col <- rainbow(length(all.data))
names(col) <- names(all.data)
for(n in names(colors)) {
    if (colors[[n]] != "")
      col[n] <- colors[[n]]
}

densities <- lapply(all.data, density)
maxy <- max(unlist(lapply(densities,function(d)max(d$y))))

out <- sub("\\.pdf$","", args$out)
out <- paste0(out, ".pdf")
title <- args$title

warning("Creating file ", out)
pdf(file = out, title = title, useDingbats = FALSE, width = 11.7, 
    height = 8.3)

scales <- as.integer(unlist(strsplit(args$scales, ",")))
stopifnot(length(scales)>0 && sum(is.na(scales))==0)

par(mfrow=c(length(scales), 1),
    mar=c(1,5,1,1))

maxx <- args$maxlen
if(maxx==0)
   maxx <- max(unlist(lapply(all.data,max)))

for(scale in scales) { 
    plot(type="n", xlab="",ylab= "",
         x=c(0, maxx), y=c(0,maxy/scale),
         xaxt="n", yaxt="n")
    s <- seq(0, maxx, 10)
    sl <- s; sl[ sl %% 50 !=0 ] <- NA
    axis(side=1,at=s, labels=sl)
    axis(side=2,labels=FALSE)
    title(ylab=sprintf("density x %.0f", scale))
    
    if (abs(scale - 1) < 1e-6) {
        title(main=title)
        legend(x="topright", legend=names(all.data),
               lty=1, col=col[names(all.data)])
    }
    for(sample in names(all.data)) {
        d <- densities[[sample]]
        d$y <-  d$y * scale
        lines(d, col=col[sample])
    }
}

dev.off()

warning("Succesfully completed\n")
print("Succesfully completed\n")
sessionInfo()
