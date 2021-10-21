\name{align}
\alias{align}

\title{
  NOT READY Align coordinates of a GRanges object. 
}

\description{
  For 'meta-gene' analyses it is often needed to align features such as
  binding sites by a set of other coordinates, such as transcription
  start sites. This function does that.
}

\usage{
  align(gr, ref = gr, width = NA, start = TRUE, ignore.strand = FALSE)
}

\arguments{
  \item{gr}{
    The features whose coordinates are to be aligned
  }

  \item{ref}{
    The reference to align to, typically the TSS of genes
  }

  \item{width}{
    If not provided, the length (and strands) of \code{gr} and \code{ref}
    have to be identical, and each range in \code{gr} is aligned by the
    start of each range in \code{ref}. This is typically used for aligning
    just genes, i.e. making their 5'-coordinates all 1, and their
    3'-coordinates equal to their own length.

    If \code{width} \emph{is} provided, it is passed to
    \code{\link{findOverlaps}} in order to match the ranges in \code{gr}
    with those in \code{ref} to do determine the amounts by which to
    shift. The ranges in \code{gr} are shifted by the starts of the
    ranges in \code{gr} that are closer than \code{width}.  }

  \item{start}{
    Also passed to \code{\link{findOverlaps}}     
  }

  \item{ignore.strand}{
    Passed to \code{\link{findOverlaps}}. For binding sites this should
    be \code{TRUE}. 
  }
}

\value{
  A new \code{GRanges} object is returned, with aligned \code{start} and
  \code{end}. The \code{strand} is set to '*'. The ordering of ranges
  may be different, and \code{width} was suppolied, the length of the
  returned object may be different, depending on what
  \code{findOverlaps()} found.  Note that \code{start} and \code{end} of
  transcription factor binding sites will often by negative.
}

\author{
   Philip Lijnzaad <plijnzaad@gmail.com>
}

\note{
  Make sure that the \code{\link{seqlevels}} between \code{gr} and
  \code{ref} are identical.

  Oh, and the code has not been tested thoroughly.
}

\seealso{
 \code{\link{GenomicRanges::findOverlaps}}
}

\examples{
\dontrun{

    library(rtracklayer)
    gff <- import("cerevisiae-nofasta-noquotes.gff")
    tfbs <- import("MacIsaac_2006_ChIP_chip_TFBSs_V64.gff3")
    plot(density(mid(ranges(align(tfbs[tfbs$binding_factor=="STE12"],
                           ref=gff[gff$type=="gene"], width=600)))))
}
}

\keyword{misc}

