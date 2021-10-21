\name{location2granges}
\alias{location2granges}

\title{
Convert a string describing a genomic location to a GRanges object.

}

\description{
  Often a human-readable string like
  "chrII:2,000,000-2,200,000"
  must be parsed to a GRanges object for internal  use. This fuction
  does that. The syntax is 

  chromosome_name | chromo:start-end | chromo:start+width | chromo:middle+-halflength
  
}

\usage{
location2granges(location, seqinfo = NULL, seqlengths = NULL)
}

\arguments{
  \item{location}{
    A string (or vector thereof) that looks

    like chromosome_name | chromo:start-end | chromo:start+width | chromo:middle+-halfwidth

    All the numbers (start, end, width, middle, halfwidth) may contain
    comma-separators. The 'start-end' specification works like
    \code{IRanges(start=s,end=e)}, whereas 'start+width' works like
    \code{IRanges(start=s,with=w)}. Lastly, 'middle+-halfwidth' is
    equivalent to
    \code{IRanges(start=middle-halfwidth,end=middle+halfwidth)} (note
    that this always yields a region of uneven width).

  }

  \item{seqlengths}{
 A named vector with sequence lengths (=chromosome sizes) per
 chromosome. If not supplied, a maximum length of 2000000000L is used
 instead, but in most downstream uses this willl break and is therefore
 warned about.
}

  \item{seqinfo}{
    Can be used instead of \code{seqlengths} argument to supply the seqlengths. (see also \code{\link{GenomeInfoDb::seqinfo}}).
  }

}

\author{
 Written by plijnzaad@gmail.com
}

\note{
  There is no syntax for the strand (yet), nor is there notation for
  kBase or MBase (yet).
}

\seealso{

  \code{\link{IRanges}}, \code{\link{GRanges}},\code{\link{Seqinfo}},
  \code{\link{rtracklayer::import}}, \code{\link{ngsutils::read.bam.region}}
}

\examples{

   location2granges(location=c("chrII",
                               "chrIII:100-200",
                               "chrIII:1,333,000+4000",
                               "chrI:100000+-3000"),
                    seqlengths=c(chrI=1000000, chrII=2000000, chrIII=3000000))

}
\keyword{misc}
