#!/bin/env Rscript

library(parseargs)

library(uuutils)

invocation <- paste(commandArgs(), collapse=" ")

## show complete information
print(R.version)
print(Sys.time())
print(invocation)

overview <- function(cat(file=stderr(),
                         "insertlen-distro.R FILENAMES\n
If a file name argument starts with somename=realfilename.insertlen, somename is used; otherwise,
the filename without extension is used"))

  
args <- parseArgs(dir=NA,               # working directory (where input and chip are taken from)
                  outdir=".",          # default: same as dir
                  out="insert-length-distro.pdf",
                  title='insert lengths',
                  .overview=overview,
                  .allow.rest=TRUE
                  )

print("Arguments:\n")
print(args)
print("\n")
libraries <- c("ngsutils", "NCIS", "rtracklayer", "IRanges", "GenomicRanges")

## load them, and show where they came from:
for (lib in libraries) { 
  library(lib, character.only=TRUE)
  cat("loaded ", lib, " from ", path.package(lib), "\n")
}

files <- args$.rest

data.name <- function(file) {
    re <- "^([^=]+)=([^=]+)$"
    if(grepl(re, file, perl=TRUE)) {
        return(unlist(strsplit(file, "=")))
    }
    name <- sub("^(.*)\\.[^.]*$", "\\1", file)
    c(name,file)
}                                       #data.name

all.data <- list()

files <-
  c("Hsf1t0-1A.insertlen", "Hsf1t0-1C.insertlen", "Hsf1t0-2A.insertlen", "Hsf1t0-2C.insertlen", "Hsf1t0-3A.insertlen", "Hsf1t0-3C.insertlen", "Hsf1t30-1A.insertlen", "Hsf1t30-1C.insertlen", "Hsf1t30-2A.insertlen", "Hsf1t30-2C.insertlen", "Hsf1t30-3A.insertlen", "Hsf1t30-3C.insertlen",
  "Hsf1t0-1B.insertlen", "Hsf1t0-1D.insertlen", "Hsf1t0-2B.insertlen", "Hsf1t0-2D.insertlen", "Hsf1t0-3B.insertlen", "Hsf1t0-3D.insertlen", "Hsf1t30-1B.insertlen", "Hsf1t30-1D.insertlen", "Hsf1t30-2B.insertlen", "Hsf1t30-2D.insertlen", "Hsf1t30-3B.insertlen", "Hsf1t30-3D.insertlen")

##files <-
##  c("Hsf1t0-1B.insertlen", "Hsf1t0-1D.insertlen", "Hsf1t0-2B.insertlen", "Hsf1t0-2D.insertlen", "Hsf1t0-3B.insertlen", "Hsf1t0-3D.insertlen", "Hsf1t30-1B.insertlen", "Hsf1t30-1D.insertlen", "Hsf1t30-2B.insertlen", "Hsf1t30-2D.insertlen", "Hsf1t30-3B.insertlen", "Hsf1t30-3D.insertlen")

## files <- files[c(4:6)]
## file <- files[1]


for(file in files) {
    name.file <- data.name(file)        #pair of name, file
    x <- read.table(name.file[2])
    all.data[[name.file[1]]] <- x[[1]]
}


densities <- lapply(all.data, density)

maxy <- max(unlist(lapply(densities,function(d)max(d$y))))

col <- rainbow(length(all.data))
names(col) <- names(all.data)

col[ grep("A$", names(col)) ] = "black"
col[ grep("B$", names(col)) ] = "blue"
col[ grep("C$", names(col)) ] = "red"
col[ grep("D$", names(col)) ] = "green4"

## out <- sub("\\.pdf$","", out)
## out <- paste0(out, ".pdf")
##pdf(file = out, title = title, useDingbats = FALSE, width = 11.7, 
##    height = 8.3, ...)

scales <- log10(c(1, 2, 5, 10, 20, 50, 100))

par(mfrow=c(length(scales), 1),
    mar=c(1,5,1,1))

## maxx <- max(unlist(lapply(all.data,max)))
maxx <- 800

for(scale in scales) { 
    plot(type="n", xlab="",ylab= "",
         x=c(0, maxx), y=c(0,maxy/(10^scale)),
         xaxt="n", yaxt="n")
    s <- seq(0, maxx, 10)
    sl <- s; sl[ sl %% 50 !=0 ] <- NA
    axis(side=1,at=s, labels=sl)
    axis(side=2,labels=FALSE)
    title(ylab=sprintf("density x %.0f", 10^scale))

    if (scale<= 1e-6) {
        title(main=title)
        legend(x="topright", legend=names(all.data),
               lty=1, col=col[names(all.data)])
    }
    for(sample in names(all.data)) {
        d <- densities[[sample]]
        d$y <-  d$y * (10^scale)
        lines(d, col=col[sample])
    }
}

dev.off()

warning("Succesfully completed\n")
print("Succesfully completed\n")
sessionInfo()
