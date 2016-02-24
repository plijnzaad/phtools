#!/usr/bin/env Rscript
## merge objects found in an .RData file into a big RData file.

library(parseArgs, quietly=TRUE)

usage <- function()warning("Usage mergeRData [--regexp '^.*$' ] --out totals.rda *.rda

Note: if the resulting *.RData contains a list whose elements are of type X (say GRanges)
rather than one big GRanges, it means that during the merging library X was not loaded.
The solution is to add more ``library(therelevantpackage)'' statements to the top of this
script.
")

args <- parseArgs(.overview=usage,
                  verbose=FALSE,
                  out=NA,
                  regexp='.*',
                  preeval='',
                  .allow.rest=TRUE)
rda.files <- args$.rest
if(is.null(rda.files)) {
    warning("No input arguments\n")
    usage()
    stop()
}

library(rtracklayer,  quietly=TRUE) # this should force loading of most of the relevant bioc. libs
  
regexps <- unlist(strsplit(args$regexp, "[,;]"))
expected <-NULL                   #based on contents of first file

final <- new.env()

check.names <- function(file, names, regexps) {
    found <- NULL
    for(re in regexps)
      found <- c(found, names[grep(re, names, perl=TRUE)])
    if(any(duplicated(found))) {
        dups <- paste(found[duplicated(found)], sep="; ")
        stop("file " ,file, "regular expressions resulted in duplicated object names, ",
             dups,"Be more specific\n")
    }
    found
}

check.sets <- function(file, expected, found) { 
    if (!setequal(expected, found))
      stop("File", file, ": expected, not found: ", setdiff(expected, found),
           "\nFound, not expected: ", setdiff(found, expected), "\n")
}


for (file in rda.files) {
    env <- new.env()
    tryCatch( load(file=file, env=env), silent=FALSE )
    contents <- ls(envir=env)
    names <- check.names(file, contents, regexps)

    if (is.null(expected)) {
        expected <- names
        for (name in names) { 
            obj <- get(name, envir=env)
            assign(name, obj, env=final)
        }
    } else {
        check.sets(file, expected, names)
        for (name in names) {
            obj <- c(get(name, envir=final), get(name, envir=env)) # the actual merging
            dummy <- length(obj)
            assign(name, obj, env=final)
        }
    }
}                                       #rda.file

save(file=args$out, list=expected, envir=final)

