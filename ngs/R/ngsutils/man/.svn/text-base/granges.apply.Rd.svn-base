\name{granges.apply}
\alias{granges.apply}

\title{

  Apply a function to a \code{numeric-Rle} over the ranges
  specified by a \code{GRanges}.

}

\description{
  In order to calculate e.g. binding statistics in a promoter, one has
  to apply a function on coverage data in one object (of type
  \code{RleList}) using the coordinates of the promoters specified by a
  \code{GRanges} object. This function does that by doing the calculation and
  adding the results as new columns to the \code{mcols} attribute of
  \code{granges}.
}

\usage{
  granges.apply(granges,
                value.list,
                seqnames = NULL,
                FUNC = standard.summary,
		ignore.strands=TRUE)
}

\arguments{
  \item{granges}{
    The \code{GRanges} object to take the coordinates from, and on which
  to base the returned object.
  }

  \item{value.list}{
    A numeric Rle-list representing the measurements that are to be
    summarized by the granges.
  }

  \item{seqnames}{
    The subset of \code{seqnames()} for which to do the summarizing. If
    \code{NULL}, the (sea)names common to \code{granges} and \code{value.list}
    are taken.
  }

  \item{FUNC}{
    A function returning an simple named numeric vector containing the
    summary data. Example functions are \code{\link{standard.summary}}
    and \code{\link{sum.positive}}. If the vector is named, those names
    will become the column names of the mcols-\code{DataFrame};
    otherwise those columns will be named "V1", "V2" etc.
    }

  \item{ignore.strand}{
    Whether to disregard the strand of the range in granges when calling
    FUNC. Currently not implemented.
  }
}

\value{ A \code{GRanges} object whose \code{\link{mcols}} attribute has
  been extended with extra column(s) representing the value(s) returned
  by \code{FUNC}.  }

\author{
 Philip Lijnzaad <plijnzaad@gmail.com>
}

\note{
  To add a column named 'auc', use a function definition such as:
  \code{ function(x){c(auc=sum(x))} }. This way, the result of
  \code{sum(x)} gets the name \code{auc}.
  
  At some point it might be nice to have a bracket-like operator for
  doing this kind of stuff.
}

\seealso{
 \code{\link{interpolate.score}}, \code{\link{standard.summary}}
}

\examples{

\dontrun{

  binding <- import.bedGraph(some.file)
  interp <- interpolate.score(granges=binding, max.dist=function(d)5*median(d))
  binding.per.bp <- coverage(interp, weight="score")

  promoters <- import.gff3(other.file)

  promoters.with.binding <- granges.apply(granges=promoters,
                                          value.list=binding.per.bp,
                                          FUNC=function(b)mean(b))
}
}
\keyword{ misc }
