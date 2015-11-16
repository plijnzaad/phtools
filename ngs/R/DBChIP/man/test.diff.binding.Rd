\name{test.diff.binding}
\alias{test.diff.binding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for differential binding
}
\description{
A null hypothesis of non-differential binding is tested at each consensus site.
}
\usage{
test.diff.binding(dat, lib.size = NULL, dispersion = NULL, common.disp = TRUE, prior.n = 10, two.sample.method = "composite.null", allowable.FC = 1.5, collapsed.quant = 0.5)
}
\arguments{
  \item{dat}{a list with the item \code{site.count}, which is a matrix of read counts at binding sites.
}
  \item{lib.size}{
a vector of library size of each ChIP sample.
}
  \item{dispersion}{
The dispersion parameter in Negative Binomial distribution. Could be a numerical value or a vector with a length of the number of consensus sites.
}
  \item{common.disp}{
logical, TRUE (use common dispersion parameter for all sites) or FALSE (use site-specific dispersion).
}
  \item{prior.n}{
a parameter regulate the degree of pooling when using site-specific dispersion (\code{common.disp}=FALSE).
}
  \item{two.sample.method}{
the method to use when comparing two condition with no replicates. The default is to test a composite null that allow certain fold change \code{allowable.FC}. Otherwise user should provide a dispersion parameter.
}
  \item{allowable.FC}{
allowable fold change when testing a composite null. Default value 1.5.
}
  \item{collapsed.quant}{
the quantile to use when testing more than two conditions without replicates. Default value is 0.5, the median.
}
}
\details{
Users are recommended to study the histogram of the $p$-values for model checking. More specifically, the $p$-values between 0.5 and 1 should be roughly uniform. 
When many replicates are available, users can also randomly split biological replicates of the same condition and perform comparisons through DBChIP using the estimated dispersion parameter to check whether the $p$-values look uniform.

}

\value{
This function return the incoming \code{dat} with new field
  \item{test.stat }{a data.frame of test statistics for testing non-differential binding at each site, include p-values and fold changes.}
}
\references{
Liang, K and Keles, S (2012) \emph{Detecting differential binding of transcription factors with ChIP-seq}, 28, 121-122.
}
\author{
Kun Liang, \email{kliang@stat.wisc.edu}
}
\seealso{
\code{\link{DBChIP}}
}