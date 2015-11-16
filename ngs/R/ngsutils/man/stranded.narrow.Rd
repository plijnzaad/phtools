\name{stranded.narrow}
\alias{stranded.narrow}

\title{
  stranded.narrow - like GenomicRanges::narrow() but
  relative to the 5'-end of a feature, rather than relative to the
  left-most coordinate.
}

\description{
  The \code{\link{flank}}
  and \code{\link{promoters}}
  function take into
  account the strand of the feature  they operate on. That is, when
  applied to a reverse strand feature, \code{flank} increases the
  \code{end} of the \code{GRanges}, rather than decreasing its \code{start}.
  This function does exactly the same for the \code{\link{narrow}} : it
  adds strand-awareness, such that the \code{start} and \code{end}
  arguments are not taken as "leftmost" and "rightmost" coordinates
  within a feature, but instead as "5'-most" and "3'-most" coordinates,
  counting from its 5'-end.
}

\usage{
  stranded.narrow(x, start = NA, end = NA, width = NA, use.names=TRUE)
}

\arguments{
  \item{x}{The GRanges object being narrowed}

  \item{start}{The start, \emph{relative to the 5' end} of the feature, to
    which to narrow. If \code{NA}, the 5'-end of the object is used.
  }

  \item{end}{The end, \emph{relative to the 5' end} of the feature, to
    which to narrow. If \code{NA}, the 3'-end of the object is used.
  }

  \item{width}{Same as for \code{narrow}. See also
    \code{\link{solveUserSEW0}}. If \code{width} is specified, exactly
    one of ( \code{start}, \code{end}) has to be specified.
  }
}

\value{Returns a narrowed GRanges object}

\author{
  Philip Lijnzaad <plijnzaad@gmail.com>
}

\note{
  Features whose strand is '*' are treated as lying on the forward
  strand; a warning is printed in this case.

  In the future there may be a GenomicRanges::narrow( ..., ignore.strand=TRUE)
   method that is strand-aware, but that defaults to the old IRanges
  behaviour.
}

\seealso{
  \code{\link{narrow}},
  \code{\link{flank}},
  \code{\link{promoters}},
  \code{\link{solveUserSEW0}}
}


\examples{

gr <- GRanges(seqnames=factor(letters[1:3]),
              ranges=IRanges(start=c(101,101,101), end=c(200,200,200)),
              strand=c('+', '-', '*'),
              seqlengths=c(a=500,b=500, c=500))

narrow(gr, start=21, end=50)

stranded.narrow(gr, start=21, end=50)

}

\keyword{misc}
