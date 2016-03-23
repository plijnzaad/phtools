\name{decommafy}
\alias{decommafy}
\title{
   Get rid of comma's (used for integer-containing strings formatted with separator comma's  for readability)
}
\description{
   \code{decommafy( "1,234,567,809" )} -> 1234567809L

}
\usage{
  decommafy(x)
}

\arguments{
  \item{x}{
    Integer-containing strings formatted with separator comma's  for readability
  }

}

\value{

}

\author{
plijnzaad@gmail.com
}

\seealso{
\code{\link{commafy}}
}

\examples{

  decommafy( c("10", "1,234,567,809", "1e6"))

}

\keyword{misc}

