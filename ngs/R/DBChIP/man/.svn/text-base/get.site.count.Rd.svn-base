\name{get.site.count}
\alias{get.site.count}
\title{
Get site count
}
\description{
Count number of reads around each binding site
}
\usage{
get.site.count(dat, subtract.input=FALSE, window.size = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a list with the following items: consensus.site, chip.list, input.list, matching.input.names, norm.factor.vec. Description of the items can be found in the return value of \code{\link{DBChIP}}.
}
  \item{subtract.input}{
logical. Whether input will be subtracted from ChIP when counting the binding reads. Default is FALSE.
}
  \item{window.size}{
the window size to count reads around each binding site. Default 250 bp.
}
}
\details{
The read count for each binding site is defined as the sum of the number of 5' ends on the positive strand within the upstream window \code{[s-w, s-1]} and the number of 5' ends on the negative strand within the downstream window \code{[s+1, s+w]}, where \code{s} is a consensus site position and \code{w} is the window size parameter.
}
\value{
This function return the incoming \code{dat} with the new component:
  \item{site.count }{a matrix of read counts for each site (row) and each ChIP library (column).}
}
\seealso{
\code{\link{DBChIP}}
}