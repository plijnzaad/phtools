\name{binding.site.list}
\alias{binding.site.list}
\docType{data}
\title{
Binding site predictions
}
\description{
Binding site predictions for emb and L1 conditions in a study of transcription factor PHA-4/FOXA in \emph{C.elegans}
}
\usage{data("PHA4")}
\format{
List of 2 elements: emb and L1. Each element is a data.frame with fields: chr, strand and weight (optional).
}
\details{
Weight represents a measure of strength of the binding, for example, the number of reads in the peak.
}
\source{
Zhong et al. (2010), \emph{Genome-wide identification of binding sites defines distinct functions for Caenorhabditis elegans PHA-4/FOXA in development and environmental response}, PLoS Genetics, 6, e1000848.
}

\examples{
data(PHA4)
names(binding.site.list)
head(binding.site.list[["emb"]])
}

\keyword{datasets}

\seealso{
\code{\link{PHA4}}
}