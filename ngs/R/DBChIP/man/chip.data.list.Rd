\name{chip.data.list}
\alias{chip.data.list}
\docType{data}
\title{
ChIP data of transcription factor PHA-4/FOXA in \emph{C.elegans}
}
\description{
This data set contains parts of ChIP-seq data of transcription factor PHA-4/FOXA in \emph{C.elegans}. 
Only ChIP data in chromosome I with position < 0.9M bp are included.
}
\usage{data("PHA4")}
\format{
List of 4 elements: "emb_rep1", "emb_rep2", "L1_rep1" and "L1_rep2". Each element is a data.frame with fields: chr, strand and pos.
}
\source{
Zhong et al. (2010), \emph{Genome-wide identification of binding sites defines distinct functions for Caenorhabditis elegans PHA-4/FOXA in development and environmental response}, PLoS Genetics, 6, e1000848.
}

\examples{
data(PHA4)
names(chip.data.list)
head(chip.data.list[["emb_rep1"]])
}

\keyword{datasets}

\seealso{
\code{\link{PHA4}}
}