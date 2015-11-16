\name{conds}
\alias{conds}
\docType{data}
\title{
conds
}
\description{
Experimental conditions of ChIP replicates in PHA4 data.
}
\usage{data("PHA4")}
\format{
  The format is:
 Factor w/ 2 levels "emb","L1": 1 1 2 2
}
\details{
The first two ChIP replicates are in embryonic (emb) condition, and the last two are in the first stage of larval development (L1) condition.
}
\source{
Zhong et al. (2010), \emph{Genome-wide identification of binding sites defines distinct functions for Caenorhabditis elegans PHA-4/FOXA in development and environmental response}, PLoS Genetics, 6, e1000848.
}
\examples{
data(PHA4)
conds
}

\keyword{datasets}

\seealso{
\code{\link{PHA4}}
}