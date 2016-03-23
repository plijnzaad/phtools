\name{commafy}
\alias{commafy}

\title{
Put comma's into a big integers to separate the powers of 1000 for readability
}
\description{
   \code{commafy(1234567809)} -> "1,234,567,809"
}

\usage{
  commafy(x, preserve.width = "common")
}

\arguments{
  \item{x}{ integers, or numerics or strings that can be converted to
    them using \code{\link{as.integer}}. This may result in truncation
  of fractions, and/or \code{NA}s for invalid syntax.
}
  \item{preserve.width}{ see argument of same name in \code{formatC}}
}


\value{  strings with all powers of 1000 separated by ',' for readability.
}

\author{
  plijnzaad@gmail.com
}

\seealso{
 \code{\link{decommafy}}
}

\examples{

  commafy( c(pi, 10, 1234567809, 1e6, 1e100))

}
\keyword{misc}
