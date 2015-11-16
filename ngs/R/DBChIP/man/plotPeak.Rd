\name{plotPeak}
\alias{plotPeak}
\title{
Plot peak
}
\description{
Plot the coverage profile of differentially bound peaks
}
\usage{
plotPeak(rept, dat, lib.size = NULL, w = 400, ext = 200, combine.rep = FALSE, cap = NULL, n.row.per.page = 6, caption = NULL)
}
\arguments{
  \item{rept}{
a data.frame as the result of \code{\link{report.peak}}.
}
  \item{dat}{
a list with the following items: chip.list, input.list, matching.input.names, norm.factor.vec. Description of the items can be found in the return value of \code{\link{DBChIP}}.
}
  \item{lib.size}{
a vector of library size of each ChIP sample. Used to scale profiles so that they are comparable.
}
  \item{w}{
half window size to plot around the binding sites.
}
  \item{ext}{
extension size for each read. Each read is extended from its 5' end by \code{ext} towards 3' end. If \code{ext} is set to 1, the minimum, per nucleotide read counts will be plotted.
}
  \item{combine.rep}{
logical, whether to combine replicates for the plot. Can be useful when there are too many replicates to plot.
}
  \item{cap}{
the maximum number of reads per nucleotide allowd to plot. 
}
  \item{n.row.per.page}{
the maximum number of tracks (rows) per page. 
}
  \item{caption}{
additional caption to appear on the title besides the location.
}
}
\seealso{
\code{\link{DBChIP}}
}