\name{read.macs2xls}
\alias{read.macs2xls}
\title{
  Read peaks called by MACS2 from the MACS2 *.xls file.
}
\description{

  The peakcalling program MACS2 (see reference) produces a number of output files, one of which
is the tab delimited file of peaks, summits, fold-enrichment and
p-values, and is (confusingly) called NAME.xls. This function reads that file,
and returns it as a \code{\link{GRanges}} object.

\code{\link{GRanges}} object.
}

\usage{
  gr <- read.macs2xls(file)
}

\arguments{
  \item{file}{
    Name of the file
  }
  \item{\dots}{
    All further arguments (typically, \code{seqlengths} and
  \code{seqinfo})  are passed to the \code{\link{GRanges}} function
  }
}

\value{
  A \code{\link{GRanges}} object. The meta data columns (i.e. \code{mcols(obj)}) are:

\code{abs_summit}, \code{pileup}, \code{minlog10.pval},
\code{fold_enrichment}, \code{minlog10.qval} and \code{name}
}

\references{
 Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137
}

\author{
<plijnzaad@gmail.com>
}

\note{
  The code was written for MACS2 version 2.1.1; if its format changes,
  this function will need changing.

  The comments at the top of xls file are of course skipped, but contain
  valuable information.
}

\examples{
###  peaks <- read.macs2xls(file)
}

\keyword{misc}% use one of  RShowDoc("KEYWORDS")

