\name{report.peak}
\alias{report.peak}
\title{
Report peaks
}
\description{
Report most significant peaks
}
\usage{
report.peak(test.res, FDR = NULL, FDR.method = "BH", n = 10, add.origin = TRUE, adaptive.threshold = c(0.05, 0.95))
}
\arguments{
  \item{test.res}{
a list with the item \code{test.stat}, which is a data.frame of test statistics for testing non-differential binding at each site, include p-values and fold changes.
}
  \item{FDR}{
the desireable false discovery rate (FDR) level.
}
\item{FDR.method}{
the method to control FDR. Default is "BH", the Benjamini and Hochberg (1995) method. Another option is "adaptive", which use the adaptive.threshold to estimate the number of true null hypotheses.
}
  \item{n}{
the top \code{n} differential peaks to return.
}
  \item{add.origin}{
logical. Whether to add peaks' origin information (how many and which conditions the consensus site are merged from) in the result.
}
  \item{adaptive.threshold}{
vector of two values between 0 and 1. The number of true null hypotheses is estimated as the number of p-values between these two values divided by the distance between these two values. Default is (0.05, 0.95). If it is set to (0.05, 1), it becomes the method recommended in Blanchard and Roquain (2009).
}
}
\details{
The default is to return the top \code{n} differential peaks. If user specify \code{FDR}, a set of peaks under the threshold will be returned instead. The FDR is computed through the classical Benjamini & Hochberg 1995 method.
}
\value{
a data.frame with with following fields
\item{chr}{chromosome.}
\item{pos}{consensus binding position.}
\item{nsig}{number of significant original binding sites that are merged into consensus site.}
\item{origin}{the origin of merged binding sites.}
\item{ori.pos}{the original positions of merged binding sites. The consensus position is a weighted average of original positions.}
\item{FC.condition name}{the fold change comparing to the first condition.}
\item{pval}{p-value for testing non-differential binding.}
\item{FDR}{the q-value.}
}
\references{
Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{J. R. Stat. Soc. B}, 57, 289-300.

Blanchard, G. and Roquain, E. (2009). Adaptive false discovery rate control under independence and dependence. \emph{Journal of Machine Learning Research}, 10, 2837-2871.
}
\seealso{
\code{\link{DBChIP}}
}