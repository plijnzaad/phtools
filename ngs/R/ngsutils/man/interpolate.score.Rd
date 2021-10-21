\name{interpolate.score}
\alias{interpolate.score}

\title{

  Interpolate the \code{\link{score()}} of a (sparse) signal down to
  base pair level.
}

\description{
  Tiling arrays have a lower density of information per kbase. 
  To compare data measured with such techniques with basepair resolution
  obtained from NGS, it can be useful to interpolate.score the sparser data.
}

\usage{

  interpolate.score(granges, seqnames = NULL,
                    wanted.dist = function(dists){median(dists)},
                    max.dist = function(dists){3*median(dists)},
                    every = 1)
}

\arguments{
  \item{granges}{
    The \code{GRanges} object whose \code{score} is to be interpolated.
    The \code{\link{seqlength}}'s of this object may not be NULL.
  }
  
  \item{seqnames}{
    The subset of seqnames for which to do the interpolation. If
    \code{NULL}, all \code{unique(seqnames(granges))} are taken.
  }

  \item{wanted.dist}{
    Passed to \code{\link{uuutils::zeroterminate.islands}}
  }

  \item{max.dist}{
    Passed to \code{\link{uuutils::zeroterminate.islands}}
  }

  \item{every}{
    Interpolate every this many base pairs.
  }
}

\details{

  To avoid interpolation over too long distances, the 'islands' in the
  data are first 'zero-terminated' on either side using
  \code{\link{uuutils::zeroterminate.islands}}. This makes sure that,
  for the simple linear interpolation done subsequently, the complete
  gap is interpolated as having the value zero.

}

\value{
  A \code{GRanges} object with, with 1-width ranges at the interpolate.scored
  points, and \code{score} the interpolated values.
}

\author{
  Philip Lijnzaad <plijnzaad@gmail.com>
}

\seealso{
  \code{\link{uuutils::zeroterminate.islands}},\code{\link{granges.apply}}
}

\note{ Don't be seduced by the luscious curves of spline interpolation.
  They are visually attractive but inappropriate as they
  overshoot wildly. Better stick to the straight and narrow, as
  in real life.
}

\examples{

  gr <- GRanges(ranges=IRanges(start=seq(1, 1000, by=100),width=1),
              score=rnorm(10),
              seqnames=factor('foo'),strand='*', seqlengths=c(foo=1000))

  interpolation <- interpolate.score(gr)

}
\keyword{misc}
