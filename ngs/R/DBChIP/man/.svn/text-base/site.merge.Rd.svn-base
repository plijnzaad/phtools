\name{site.merge}
\alias{site.merge}
\title{
Cluster close-by sites into consensus sites
}
\description{
Cluster close-by sites into consensus sites through agglomerative (bottom-up) hierarchical clustering.
}
\usage{
site.merge(bs.list, in.distance = 100, out.distance = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bs.list}{
internal format of binding site. Result of \code{\link{read.binding.site.list}}.
}
  \item{in.distance}{
within cluster distance. If the distance between centroids of two clusters are smaller than this value, the clusters will be merged into one. Default value 100 bp.
}
  \item{out.distance}{
out of cluster distance. If the distance between centroids of two clusters are larger than this value, they are considered different binding sites. Default value 250 bp.
}
}
\details{
We group predicted locations from multiple conditions into clusters of close-by locations 
by using agglomerative (bottom-up) hierarchical clustering with centroid linkage. 

If the distance between centroids of two clusters are between \code{in.distance} and \code{out.distance}, the cluster with higher weight will be kept.
}
\value{
A list where each element represents a chromosome. Within each chromosome, it is a data.frame with following fields
\item{pos}{consensus binding position.}
\item{nsig}{number of significant original binding sites that are merged into consensus site.}
\item{origin}{the origin of merged binding sites.}
\item{ori.pos}{the original positions of merged binding sites. The consensus position is a weighted average of original positions.}
}

\seealso{
\code{\link{DBChIP}}
}
