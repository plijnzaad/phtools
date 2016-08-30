library(Biobase)

### library(Rsamtools)

library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(ShortRead)

if(packageVersion("GenomicRanges") < "1.10.7") {
  warning("**** Update the bloody GRanges package to at least 1.10.7 ***")
  mcols <- function(grange)DataFrame(grange@elementMetadata@listData)
  ### needed by the import 
}

stranded.narrow <- function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
  ## as IRanges::narrow(), but start and end are interpreted relative to
  ## the proper 5'-end of the GRange, i.e. it does take into account the
  ## strand of the feature (which is unknown to IRanges objects).  In
  ## the future there may be a GenomicRanges::narrow( ..., ignore.strand=TRUE)
  ## that knows strands but defaults to the old IRanges behaviour

  stopifnot(is(x, "GRanges"))
  w <- width(x)
  lg <- length(x)
  strands <- strand(x)
  if(any(strands=='*'))
    warning("Some elements have no strand; they are treated as positive-stranded")
  iranges <- solveUserSEW(refwidths=width(x),
                          start=start,end=end,width=width)
  li <- length(iranges)
  if( li != lg ) {   #can't deal with anything more complicated now
    stopifnot(li == 1)
    starts <- rep(start,lg)
    ends <- rep(end,lg)
  } else { 
    starts <- start(iranges)
    ends <- end(iranges)
  }
  rev <- as.logical(strands == '-' )    # as.vector should also work?!
  orig.starts <- starts
  starts[rev] <-   w[rev] - ends[ rev ] + 1
  ends[rev] <-   w[rev] - orig.starts[ rev ] + 1
  narrow(x, starts, ends, use.names=use.names)
}                                       #stranded.narrow

read.chromSizes <- function(file) {
  x <- read.table(file)
  sizes <- x[[2]]
  names(sizes) <- x[[1]]
  sizes
}

read.complete.bigwig <- function(file, chromSizes) {
  chromSizes <- read.chromSizes(chromSizes)

  args <- list()
  for(chr in names(chromSizes))
    args <- c(args, IRanges(start=1, end=chromSizes[chr]))
  names(args) <- names(chromSizes)
  ## BigWigSelection(IRanges::RangesList(I = IRanges(start=1, end=230218)))
  genome <- BigWigSelection(IRanges::RangesList(args)) #all

  all <- import.bw(con=file, selection=genome)
}                                       #read.complete.bigwig

## .rename.chromos <- function(s)gsub("mt", "Mito", gsub("chr", "", s))
.rename.chromos <- function(s)s

read.sgd.genes <- function(...)
  stop("obsolete, read.sgd.genes() now called read.sgd.features()")

read.sgd.features <- function(file,
                           feature.type=c("gene"),
                           ### check 
                           skip.dubious=TRUE) {
  stopifnot(is.null(feature.type) || is.character(feature.type))
  if(is.null(feature.type) || feature.type == "" || feature.type == "all" )
    feature.type <- NULL                #meaning 'all'
  else
    feature.type <- union(feature.type, "chromosome") # chromo needed to get full length
  
  stopifnot(packageVersion("GenomicRanges") >= "1.6.7") # wild guess

  stopifnot(file.exists(file))

  cat("Reading ", file, "\n")

  features <- import.gff3(file, feature.type=feature.type,
                         genome="S.cerevisiae")

  ## keep the chromosome lengths:
  chromos <- features[ features$type=='chromosome',"ID"]
  chr.lengths <- width(chromos)
  names(chr.lengths) <- .rename.chromos(seqlevels(chromos))
  rm(chromos)

  features <- features[ mcols(features)$type != "chromosome" ]

  if(skip.dubious) {
    skip <- mcols(features)$orf_classification == "Dubious" # 
    skip[ is.na(skip) ] <- FALSE                            # orf_classification can be NA
    features <- features[! skip ]
  }
  names <- mcols(features)$ID
  stopifnot(all(!is.na(names) & nchar(names)>2))
  names(features) <- names

  seqlengths(features)[names(chr.lengths)] <- chr.lengths #they get lost ...
  cat("Read ", length(features), " features\n")
  features
}                                       #read.sgd.features

read.sgr <- function(file) {
    ## import an sgr file as a GRanges
    with(read.table(file, sep="\t",
                    col.names=c("seqnames", "start", "score")),
         GRanges(seqnames=seqnames,
                 ranges=IRanges(start=start, width=1),
                 strand='*',
                 score=score))
}                                       #read.sgr


center <- function(reads, chromsizes, shift) {
### shifts all reads by shift in their 3'-direction
  stop("ngsutils::center: only meaningfull if reads are reduced to 1 bp reads--fix this. PL")
  stopifnot(class(reads)=="GRanges")
  stopifnot(is.numeric(chromsizes))
  stopifnot(is.numeric(shift))

  ranges(reads[(strand(reads)=='-')]) <- shift(ranges(reads[(strand(reads)=='-')]), -shift)
  { falling.off <- strand(reads)=='-' & end(reads)<1 
    reads <- reads[! falling.off]
    cut <- strand(reads)=='-' & start(reads)<1 
    start(reads[cut]) <- 1
  }

  ranges(reads[(strand(reads)=='+')]) <- shift(ranges(reads[(strand(reads)=='+')]), shift)
  { falling.off <- rep(FALSE, length(reads))
    for(chr in names(chromsizes)) { 
      max <- chromsizes[chr]
      falling.off <-  (falling.off | seqnames(reads)== chr & strand(reads)=='+' & start(reads)> max)
    }
    reads <- reads[!falling.off]
  }
  
  for(chr in names(chromsizes) ) { 
    max <- chromsizes[chr]
    idx <- seqnames(reads) == chr & strand(reads)=='+' & end(reads) > max
    if(any(idx))
      end(reads[idx]) <- max
  }
  reads
}                                       #center

get.coverage <- function(file, chromsizes, shift=NULL) {
    ## reads .bed file, maybe centers by shift bp in 3'-direction, and
    ## returns coverage. 
    warning("**** NOTE: does not yet deal with paired end data! ****\n")
    cat("reading ", file, "\n")
    reads <-import(file, format="bed")
    if (!is.null(shift)) { 
        cat("centering ", file, "\n")
        reads <- center(reads=reads, chromsizes=chromsizes, shift=shift)
    }
    cat("calc. coverage ", file, "\n")
    coverage(reads)
}                                       #get.centered.coverage

igb.format <- function(grange) {
  ### utility to format the coordinates of a GRanges object
  ### for copy-pasting into the coordinates field of IGB or igv
  stopifnot(is(grange, "GRanges"))
  chr <- as.character(seqnames(grange))
  start <- commafy(as.vector(start(grange)))
  end <- commafy(as.vector(end(grange)))
  paste(" ", gsub(" ","", sprintf("%s:%s-%s", chr, start,end))," ")
}

.codname <- function(file)gsub("\\.[^.]*$", "", basename(file))
  ### unlist(strsplit(basename(file), "\\."))[1]

read.cod <- function(file, prefix= .codname(file), seqlengths=NULL,
                     score.col="max_log2FC"
                     ## Which column to use as the score of the feature.
                     ## Other possibilities are "FDR" , "bound_width" ,
                     ## "maxT", , "minuslog10_minPoisP" or any of the
                     ## "normalized DNA fragment counts" columns (which
                     ## may have gotten weird names during the import).
                     ) {
### CisGenome produces a big output table called something.cod. This function
### reads it in as a GRanges object. The prefix arg is prepended to the peak rank
### to create the seqname. For  information on the meaning of all the columns,
### see http://www.biostat.jhsph.edu/~hji/cisgenome/index_files/outputformat_seqpeak.txt
### 
    
    x <- read.table(file = file, sep = "\t", as.is = T, quote = "", header = T, 
                    comment.char = "", row.names = 1)
    nrow <- nrow(x)

    if (! score.col %in% names(x))
        stop("Column ", score.col, " not found, cannot use to assign scores")

    warning("Usings column ", score.col, " as the score for these peaks\n")
    
    arg.list <- list(seqnames=Rle(x$chromosome),
                     ranges=IRanges(start=x$start,
                         end=x$end,
                         names=paste(prefix,rownames(x),sep="")),
                     strand=Rle('*', nrow),
                     seqlengths=seqlengths,
                     score=x[[score.col]])

    ## add remaining columns as mcols
    ##  ignore <- c("chromosome", "start", "end", "strand", "peak_length", score.col)
    ## no, keep score.col for clarity
    ignore <- c("chromosome", "start", "end", "strand", "peak_length")
    mcol.names <- setdiff(names(x), ignore)
    for(mcol.name in mcol.names) { 
        col <- list(x[[ mcol.name ]])
        names(col) <- mcol.name
        arg.list <- c(arg.list, col)
    }

    do.call(GRanges, arg.list)
}                                       #read.cod

subsetByBin <- function(granges, bin, ignore.strand=FALSE) {
  ## give the GRanges subset (boolean vector) that corresponds to given bin.
  ## typical usage: 
  ## bins <- IRanges(start=seq(from=1, to=len,by=bin.width),
  ##                 width=bin.width)
  ## means <- sapply(1:length(bins), function(i)mean(score(chr[ subsetByBin(chr, bins[i] ) ])))
  ## NOTE: see also binnedAverage() described in the tileGenome man page!

    
  stopifnot(is(granges,"GRanges"))
  stopifnot(is(bin,"IRanges"))
  subset <- as.logical(findOverlaps(ranges(granges), bin, select='first'))
  ### NOTE: this ignores the fact that the width of the is not always 1,
  ### but that makes it much easier
  subset[is.na(subset)] <- FALSE ## findOverlaps returns 1 and NA, not TRUE and FALSE
  subset
}                                       #subsetByBin

binwise.aggregate <- function(granges, bins, FUN=mean, empty.bin.value=0) {
  stopifnot(is(granges, "GRanges"))
  stopifnot(is(bins, "IRanges"))
  x <- sapply(1:length(bins), function(i)FUN(score(granges[ subsetByBin(granges, bins[i] ) ])))
  sane <- is.finite(x) & (!is.nan(x)) & !(is.na(x))
  x[!sane] <- empty.bin.value
  x
}                                       #binwise.aggregate

intergenic.features <- function(features) {
  stopifnot(is(features, "GRanges"))
  ## returns convergent intergenic regiongs as a set of features (labeld by
  ## name of ORF to the left of it
  grl <- GRangesList()                  #ordinary list won't work

  for(chr in unique(seqnames(features))) { 

    genes <- features[ seqnames(features) == chr ]
    p <- precede(genes, genes, ignore.strand=TRUE)
### stopifnot(is.na(p[length(p)]) && sum(is.na(p))==1)
    ignore <- is.na(p)
    p <- p[ !ignore ]
    follow <- genes[p]
    genes <- genes[ !ignore]
    stopifnot(sum(end(genes) >= start(follow))==0)

    subsets <- list(tandem=strand(genes)== strand(follow),
                    convergent=strand(genes)=="+" & strand(follow)=="-",
                    divergent=strand(genes)=="-" & strand(follow)=="+")

    for(name in names(subsets))
      if(sum(subsets[[name]])==0)       # e.g. if '*' or for mitochondria ...
        subsets[[name]] <- NULL

    intergenic <- IRanges(start=end(genes)+1, end=start(follow)-1,
                          names=paste(sep="","i",names(genes)))
    sets <- GRangesList()                  #ordinary list won't work
    for(name in names(subsets)) {
      g <- GRanges(seqnames=chr, ranges=intergenic[ subsets[[name]] ], strand='*')
      g$type <- sprintf("intergenic_%s", name)
      sets[[name]] <- g
    }
    grl[[chr]] <- unlist(sets)
  }                                     #for chr
  unlist(grl)
}                                       # intergenic.features

zero.fill <- function(track) { 
  ## make a GRange covering the whole genome, including all the zeros.
  ## For now, assume one chromosome and strand:
  chr <- unique(seqnames(track))
  stopifnot(length(chr)==1)
  strand <- unique(strand(track))
  stopifnot(length(strand)==1)
  g <- gaps(track)
  g <- g[ seqnames(g)==chr & strand(g)==strand ]
##  values(g)$score <- Rle(0)                  #score?
  score(g) <- Rle(0)
  sort(c(g, track))
}                                       #zero.fill

### straight from the IRanges vignette (better use Gviz)
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, add=FALSE, ...) {
    height <- 1
    if (is(xlim, "Ranges"))
        xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    if (!add)  {
        plot.new()
        plot.window(xlim, c(0, max(bins)*(height + sep)))
        title(main)
        axis(1)
    }
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
}                                       # plotRanges

interpolate.score <- function(granges,         # a GRanges
                              seqnames=NULL,
                              wanted.dist = function(dists){median(dists)},
                              max.dist = function(dists){3*median(dists)},
                              every=1) {                    #how/where to interpolate
### given a signal from a tiling array, interpolate it such that we have
### a signal for each bp, i.e. pretty much like a chipseq track
    library(uuutils)
    stopifnot(is(granges, "GRanges"))
    l <- seqlengths(granges)
    if( is.null(l) ||  any(is.na(l)) )
        stop("interpolate.score: seqlengths() of granges argument must have been set")
    stopifnot(is.numeric(score(granges)))
    if(is.null(seqnames))
        seqnames <- unique(as.character(seqnames(granges)))
    new.granges <- list()
    for (seq in seqnames ) {
        s <- granges[ seqnames(granges)==seq ] # s for 'signal'

        ## actual interpolation. First order the separate GRanges segments by their midpoints:
        x <- floor((start(s) + end(s))/2)   # in case width!=1
        o <- order(x)
        x <- x[o]
        y <- (score(s))[o]
        ## surround the gaps with 0's:
        zeroterm <- uuutils::zeroterminate.islands(x,y,
                                                   wanted.dist=wanted.dist,
                                                   max.dist=max.dist) 
        len <- seqlengths(s)[seq]
        n.out <- floor(len/every)
        interp <- approx(zeroterm$x, zeroterm$y,
                         xout=floor(seq(1, len, length.out=n.out)),
                         yleft=0,yright=0)
        gr <- GRanges(seqnames=seq,
                      ranges=IRanges(interp$x, width=1), # singe bp irange at the midpoint
                      strand='*',
                      seqlengths=seqlengths(granges)[seq],
                      score=interp$y
                      ## seqlevels=names(seqlengths(granges)) ### BROKEN?!
                      )
        seqlevels(gr) <- names(seqlengths(granges))
        new.granges <- c(new.granges, gr)
    }   #seq
    all <- do.call("c", args=new.granges)
    seqlevels(all) <-  as.character(unique(seqnames(all)))
    ## (this is to avoid the
    ## "Error in .Call2("Rle_constructor", values, lengths, check, 0L, PACKAGE = "IRanges") : 
    ##     'length(lengths)' != 'length(values)'. (Fixed in newer versions?)
    all
}                                       #interpolate.score

standard.summary <- function(val){c(sum=sum(val, na.rm=TRUE),
                                     min=min(val, na.rm=TRUE),
                                     mean=mean(val, na.rm=TRUE),
                                     median=median(val, na.rm=TRUE),
                                     max=max(val, na.rm=TRUE),
                                     sd=sd(val, na.rm=TRUE),
                                     IQR=IQR(val, na.rm=TRUE))}

sum.positive <- function(val){c(sum.pos=sum(val[val>0], na.rm=TRUE))}


granges.apply <- function(granges,
                          value.list,
                          seqnames=NULL,
                          FUNC=standard.summary,
                          ignore.strand=TRUE) {
    ## to granges, add columns that summarize (using the specified
    ## FUNC) the values found in value.list, over the ranges
    ## of the granges. The value.list is an RleList that typically came
    ## from calling coverage().

    stopifnot(is(granges, "GRanges"))
    stopifnot(is(value.list, "List")|| is(value.list, "list"))
    ## from e.g. coverage()
    ## how can we check contents that are numeric??
    granges.names <- unique(as.character(seqnames(granges)))
    if(is.null(seqnames))
        seqnames <- intersect(granges.names, names(value.list))
    if(is.null(seqnames))
        stop("granges.apply: no seqnames. Were they set for the granges argument, and was there at least some overlap with names(value.list)?")
    granges.names.notfound <- setdiff(granges.names,names(value.list))
    values.names.notfound <- setdiff(names(value.list), granges.names)
    if(length(granges.names.notfound)>0)
        warning(paste("granges.apply: coordinates, but no values found for following seqnames:", granges.names.notfound, collapse=" "))
    if(length(values.names.notfound)>0)
        warning(paste("granges.apply: values, but no coordinates found for following seqnames:", values.names.notfound, collapse=" "))

    if(!ignore.strand)
        stop("Not yet implemented")
    new.granges <- list()
    for(seq in seqnames) {
        feats <- granges[ as.character(seqnames(granges))==seq ]
        vals <- as.numeric((value.list[[seq]])[ranges(feats)])
        ## Note: this is all concatenated, disentangle:
        w <- width(feats)
        start.end <- cbind( cumsum(c(1, w[-length(w)])), cumsum(w))
        vals.per.feat <- apply(start.end,1, function(pair)vals[pair[1]:pair[2]])
        ## NOTE: result may be a list or matrix, act accordingly:
        if (is.list(vals.per.feat)) {
            summary <- lapply(vals.per.feat,FUNC) # list of named vectors 
            colnames <- names(summary[[1]]) ## "sum", "min",etc
            if(is.null(colnames))
                colnames <- paste0("V", 1:length(summary[[1]]))
            summary <- matrix(unlist(summary),byrow=TRUE,
                              nrow=length(feats),
                              ncol=length(colnames),
                              dimnames=list(names(feats),colnames))
            ## rows=features, columns = result of FUNC, named 
        } else {  # not a list, but simplified to a plain matrix:
            summary <- t(apply(vals.per.feat,2, FUNC))
            rownames(summary) <-  names(feats)
        }
        new <- c(mcols(feats), DataFrame(summary))
        stopifnot(is(new, "DataFrame")) # double check just in case ...
        mcols(feats) <- new
        new.granges <- c(new.granges, feats)
    }
    do.call("c", new.granges)
}                                       #granges.apply

read.bam.region <- function(file, region=NULL) {
    if(is.null(region) || length(region)==0 || is.na(region) || region=="")
      return(as(import(file, format="bam"), "GRanges"))

    ### parts <- unlist(strsplit(gsub("[,_]","", region), "[-:]"))
    parts <- sapply(region, function(r)unlist(strsplit(gsub("[,_]","", r), "[-:]")))
    if(is.list(parts))
      stop("Whens specifying several regions they must all have the same number of components")
    if(!is.matrix(parts))               #one or more full chromosomes
      parts <- rbind(seq=parts, start=1, end=536870912L)
      ##(weird maximum, just twice human chrI, not enough for Paris Japonica ...)

    if (nrow(parts)==2)                 # chrX:1000 means: from 1000 till end ...
      ## stop("Invalid region syntax: ", region, ", expected seqName:start-end\n")
      parts <- rbind(parts, end=536870912L)

    seq <- parts[1,]
    dim <- dim(parts[2:3,,drop=FALSE])
    parts <- as.integer(parts[2:3,])
    dim(parts) <- dim                # ,drop=FALSE does not help, as.integer still drops it
    rownames(parts) <- c('start', 'end')

    if(any(is.na(parts)))
      stop("Invalid region syntax: ", region, ", expected seqName:start-end\n")

    param <- ScanBamParam(flag=scanBamFlag(),
                          simpleCigar=TRUE,
                          reverseComplement=FALSE,
                          which=GRanges(seq,IRanges(parts['start',], parts['end',])),
                          what=c('rname', 'pos', 'strand', 'qname'))

    reads <- readGAlignments(file=file, param=param)
    GRanges(seqnames=seqnames(reads),
            ranges=ranges(reads),
            strand=strand(reads),
            seqinfo=seqinfo(reads))
}                                       #read.bam.region


align <- function(gr, ref=gr, width=NA, start=TRUE, ignore.strand=TRUE) {
### first attemp to do proper meta-gene plots
    warning("not yet ready")
    stopifnot(is(gr, "GRanges"))
    stopifnot(is(ref, "GRanges"))
    stopifnot(identical(seqlevels(gr), seqlevels(ref)))
    
    if(is.na(width)) {
        ## simple shifting, e.g. of promoters.  First some checks to see
        ## if the length and ordering is identical.  Very tricky this,
        ## will miss many cases.
        if (length(gr)!=length(ref))
            stop("GRange to be aligned not same length as that of reference ")
        if(!all(seqnames(gr)==seqnames(ref)))
            stop("GRange to be aligned not same names as that of reference ")
        if(!all(   strand( gr[ strand(gr)!= '*'] )
                == strand(ref[strand(ref)!='*' ])))
            stop("Mismatch in strandedness for features to be aligned")
        if(! isTRUE(identical(names(gr), names(ref))) ) { 
            if(is.null(names(gr)))
                warning("Neither gr nor ref have names(), cannot check if ordering is correct")
            else
                warning("Mismatch in names of gr and ref argument, may have been reordered?")
        }

        new <- gr
        values(new) <- DataFrame(values(gr), #keep the old coords too
                                 orig.start=start(gr),
                                 orig.end=end(gr),
                                 orig.strand=strand(gr))
        ## do separately for fwd and rev
        fwd <- strand(new) %in% c("%", "+")
        ranges(new[ fwd ]) <- shift(use.names=TRUE,
                                    ranges(new[fwd]), -start(ref[fwd])+1)
        # c'est tout for fwd, now rev:
        rev <- strand(new) == "-"
        new.rev.ranges <- shift(use.names=TRUE,
                                ranges(new[rev]), -end(ref[rev]))
        start(new[rev]) <- -end(new.rev.ranges) +1
        end(new[rev]) <-  -start(new.rev.ranges) +1

        strand(new) <- Rle('*')

        return(new)
    }                                   #is.na(width)
    ## else: first match regions with the relevant ref parts, then shift
    fl <- flank(ref, width=width, start=start)
    values(fl) <- DataFrame(orig.start=start(ref),
                            orig.end=end(ref),
                            orig.strand=strand(ref))
    ref.fwd <- fl[strand(fl)  %in% c('*','+') ]
    o <- findOverlaps(gr, ref.fwd, ignore.strand=ignore.strand)

    new.fwd <- gr[ queryHits(o) ]
    values(new.fwd) <- DataFrame(values(new.fwd),
                                 orig.start=start(new.fwd),
                                 orig.end=end(new.fwd),
                                 orig.strand=strand(new.fwd))
    ranges(new.fwd) <- shift(use.names=TRUE,
                             ranges(new.fwd),
                             -ref.fwd[subjectHits(o)]$orig.start+1)
    rm(ref.fwd)
    strand(new.fwd) <- Rle('*')
    # c'est tout for the fwd, now rev
    ref.rev <- fl[strand(fl)  == '-' ]
    o <- findOverlaps(gr, ref.rev, ignore.strand=ignore.strand)

    new.rev <- gr[queryHits(o)]
    values(new.rev) <- DataFrame(values(new.rev),
                                 orig.start=start(new.rev),
                                 orig.end=end(new.rev),
                                 orig.strand=strand(new.rev))
    new.rev.ranges <- shift(use.names=TRUE,
                            ranges(new.rev),
                            -ref.rev[subjectHits(o)]$orig.end)
    start(new.rev) <-  -end(new.rev.ranges)+1
    end(new.rev) <- -start(new.rev.ranges)
    both <- c(new.fwd,new.rev)
    strand(both) <- Rle('*')
    both
}                                       #align


decommafy <- function(x){
  res <- as.integer(gsub(",", "",x))
  dim(res) <- dim(x)
  res
}

commafy <- function(x, preserve.width="common") { 
  formatC(as.integer(x), format="d", big.mark=",", preserve.width=preserve.width)
}

location2granges <- function(location, seqinfo=NULL, seqlengths=NULL) {
    ### syntax: location=loc;loc;loc
    ### loc= completechromo | chromo:start-end | chromo:start+length | chromo:start+-halflength
    stopifnot(is.character(location))

    .complete.re <- function(s)paste0("^", s, "$")
    .combine <- function(one, other, loc) {
        if (length(one) ==0 && length(other)==0)
          stop("Location ", loc, " has wrong syntax")
        if (length(one) >0 && length(other)>0)
          stop("Programmming error for (error in regexp?) ", loc, ": ", paste(collapse="\n", c(strwrap(one), strwrap(other))))
        c(one, other)
    }

    .find.coords <- function(v) {
        chr <- v[2]
        chrlen <- NA

        .check.end <- function(end, max=NULL) {
            if (is.null(max) || is.na(max)) max <- 2000000000L
            if (is.na(end)) {
                warning("Length for chromosome ", chr, " unknown, using ", max, " instead, may fail")
                return(max)
            }
            if (end>max) {
                warning("Requested end larger than chromosome length ",max,"; using latter instead.\n")
                return(max)
            }
            end
        }

        if(!is.null(seqlengths) )
          chrlen <- seqlengths[chr]

        if(length(v)==2)
          return(data.frame(chr=chr, start=1L,  end= .check.end(end=chrlen)))
            
        if(length(v)!=5)
          stop("Should not happen, error in regexp?", v[1])

            op <- v[4]
        se <- as.integer(decommafy(v[c(3,5)]))
        if(any(is.na(se)))
          stop("wrong syntax for coordinates (or perhaps too large): ", v[1] )
        if (op=='-')
          return(data.frame(chr=chr, start=se[1],  end=.check.end(se[2], max=chrlen)))
        if (op=='+')
          return(data.frame(chr=chr, start=se[1],  end= .check.end(se[1]+se[2]-1, max=chrlen)))
        if (op=='+-')
          return(data.frame(chr=chr, start=max(1L, se[1]-se[2]),  end=.check.end(se[1]+se[2], max=chrlen)))
        stop("Should not happen, error in regexp?", paste(collapse="\n", strwrap(v)))
    }                                   #.find.coords

    .check.seqnames <- function(wanted, seqnames) {
        d <- setdiff(wanted, seqnames)
        if(length(d)>0)
          stop("These chromosome(s) are not part of this genome: ", paste(d, collapse=", "), "\n")
    }                                   #.check.seqnames

    if(!is.null(seqinfo) && is.null(seqlengths))
      seqlengths <- seqlengths(seqinfo)
    
    chr.re <- '([a-z0-9]+)'
    frag.re <- '([0-9,]+)([-+]{1,2})([0-9,]+)'
    frag.re <- .complete.re(paste0(chr.re, ":", frag.re))
    chr.re <-  .complete.re(chr.re)

    chrs <- regmatches(location, regexec(chr.re, location, ignore.case=TRUE))
    frags <- regmatches(location, regexec(frag.re, location, ignore.case=TRUE))
    all <- mapply(.combine, chrs, frags, location,
                  SIMPLIFY=FALSE)# list[[1:n]] of c(inputstring, chr, location, operator, number)

    l <- lapply(all, .find.coords)
    chr <- sapply(l, function(elt)elt$chr)
    start <- sapply(l, function(elt)elt$start)
    end <- sapply(l, function(elt)elt$end)
    if(!is.null(seqinfo))
      .check.seqnames(chr,seqnames(seqinfo))
    if(!is.null(seqlengths))
      .check.seqnames(chr,names(seqlengths))
    GRanges(seqnames=chr,ranges=IRanges(start,end),
                 seqinfo=seqinfo, seqlengths=seqlengths)
}                                       #location2granges

read.macs2xls <- function(file, seqlengths=NULL, ...) {
    peaks <- read.table(file, header=TRUE)
    mcols <- peaks[ , setdiff(names(peaks), c('chr', 'start', 'end', 'length'))]
    colnames(mcols)[c(3,5)] <- paste0('minlog10.', c("p","q"), "val")
    with(peaks,GRanges(ranges=IRanges(start=start, end=end), seqnames=chr, strand='*', mcols, ...))
}

..53prime <- function(gr, side) {
    stopifnot(is(gr, "GRanges"))
    if(any(strand(gr) == "*"))
      stop("got unstranded feature")
    start(resize(gr,width=1,fix=ifelse(side==5,'start','end')))
}

#' Get the 5' (or 3') coordinates from a GRanges object.
#' @param gr The GRanges on which to operate
#' @seealso \code{start, end, width, strand}
fiveprime <- function(gr)..53prime(gr, side=5)
#' @rdname fiveprime
threeprime <- function(gr)..53prime(gr, side=3)

### Saved here too (also in /home/gen/philip/git/phtools/ngs/R/ngsutils/R/utilities.R, git version 29aug2016, but
## may disappear as it is actually close to  reduce)
#' Fuse two sets of GRanges into larger granges
#' 
#' When comparing sets of features such as binding peaks, it is
#' sometimes necessary to fuse their GRanges into larger GRanges.
#' This is made difficult by the fact that often such ranges are clustered.
#
#' This function constructs 'islands' of connected components and fuses them
#' one GRange element for all elements of the cluster. See example
#' 
#' @param query,subject,... as for \code{findOverlaps}
#' @param keep.singletons if FALSE, do not return GRange elements that did not overlap with anything
#' @return a GRanges object. \code{n.orig.members} contains the number of 
#' original GRanges elements fused into the current one. The
#' \code{n.queryhits} shows how many came from \code{query},
#' \code{n.subjecthits} shows how many came from \code{subject}. If
#' \code{keep.singletons} was TRUE, their sum is always
#' \code{n.orig.members}
#' @examples
#' gr1 <- GRanges(ranges=IRanges(start=c(1,21,31),width=c(5,5,12)), strand='*',seqnames='X', mcols=DataFrame(ID=letters[1:3]))
#' gr2 <- GRanges(ranges=IRanges(start=c(11, 19, 41),width=c(5,14,12)),strand='*',seqnames='X', mcols=DataFrame(ID=LETTERS[1:3]))
#' fuseOverlaps(gr1,gr2)
#' @note In this simple case, \code{fuseOverlaps(gr1,gr2)} is identical to \code{reduce(gr1,gr2)}
fuseOverlaps <- function(query, subject, keep.singletons=TRUE, ...) {
    library(igraph)
    overlaps <- findOverlaps(query, subject, ...)

    singles <- c(1:length(query), -(1:length(subject)))
    if(!keep.singletons) singles <- NULL
    m <- as.matrix(overlaps)
    m[,2] <- -m[,2]                     # weird bipartite graph, encode subjectHits as negative numbers

    g <- graph.empty(directed=FALSE) +  vertices(unique(c(singles, as.vector(m))))
    for(i in 1:nrow(m))
      g <- g +  edge( V(g)[[ as.character(m[i,1]) ]], V(g)[[ as.character(m[i,2]) ]])
    clusters <- components(g)$membership #named vector: names are hits, values are cluster id
    d <- data.frame(hit=names(clusters), cluster=as.character(unname(clusters)))
    ## convert to list so that clusters[[clustername]] contains all its members
    clusters <- unstack(d, form=hit~cluster)

    .get.range <- function(idx) if(idx > 0 ) query[idx] else subject[ - idx ] 

    gr <- GRanges()                   #accumulator
    for(clus in clusters) {
        fuse <- do.call(range, args=sapply(as.integer(clus), .get.range))
        h <- as.integer(clus)
        values(fuse)$n.orig.members <- length(h)
        q <- h[h>0]; s <- h[h<0]
        values(fuse)$n.queryhits <- length(q)
        values(fuse)$n.subjecthits <- length(s)
        gr <- c(gr, fuse)
    }
    gr
}                                       # fuseOverlaps
