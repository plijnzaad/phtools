\name{standard.summary}
\alias{sum.positive}
\alias{standard.summary}

\title{
  Functions to be used as the \code{summary.function}-argument of
  \code{\link{summarize.by.granges}}
}

\description{

  The \code{\link{summarize.by.granges}} function takes a functional
  argument which is called to summarize the data in \code{values.list}
  argument of \code{summarize.by.granges} of that are bounded by the
  coordinates found in the \code{granges} argument of the same.

}

\usage{
  summarize.by.granges(gr, my.values, summary.function=standard.summary)
  summarize.by.granges(gr, my.values, summary.function=sum.positive)
}

\arguments{
  \item{val}{The values to be summarized}
}

\value{
  A named vector is returned. For \code{sum.positive}, it contains just
  the sum of the positive values; for \code{standard.summary}, it
  contains the  sum, minimum, mean, maximum, standard deviation and
  interquartile range.
}

\note{
  All \code{NA} values are filtered out first.
}

\author{
   Philip Lijnzaad <plijnzaad@gmail.com>
}

\seealso{
  \code{\link{summarize.by.granges}}
}


\examples{

  standard.summary(rcauchy(100))

}

\keyword{misc}
