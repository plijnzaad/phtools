\name{PHA4}
\docType{data}
\alias{PHA4}
\title{
ChIP-seq data of transcription factor PHA-4/FOXA in \emph{C.elegans}
}
\description{
This data set contains parts of ChIP-seq data of transcription factor PHA-4/FOXA in \emph{C.elegans}. 
ChIP and control data and identified binding sites in chromosome I with position < 0.9M bp are included.
}
\usage{
data("PHA4")
}
\format{
PHA4 has three elements:
\itemize{
\item \code{binding.site.list} a list of binding sites for embryonic and L1 conditions
\item \code{chip.data.list} a list of ChIP reads data for embryonic and L1 conditions
\item \code{input.data.list} a list of control reads data for embryonic and L1 conditions
\item \code{conds} a vector of conditions of ChIP libraries.
}
}
\references{
Zhong et al. (2010), \emph{Genome-wide identification of binding sites defines distinct functions for Caenorhabditis elegans PHA-4/FOXA in development and environmental response}, PLoS Genetics, 6, e1000848.
}
\seealso{
\code{\link{DBChIP}}
}