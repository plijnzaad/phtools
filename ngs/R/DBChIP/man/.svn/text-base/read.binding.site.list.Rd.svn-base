\name{read.binding.site.list}
\alias{read.binding.site.list}
\title{
read.binding.site.list
}
\description{
Process the list of binding sites into an internal format
}
\usage{
read.binding.site.list(binding.site.list)
}
\arguments{
  \item{binding.site.list}{
a list of data.frames. Each data.frame corresponds to one condition in comparison 
and has three fields, chr, pos and weight, to indicate the binding location and strength.
}
}
\details{
The binding site for each condition is split into a list where binding sites are grouped according to their chromosomes.
}
\value{
A list of conditions, each of which is a list of binding sites in a certain chromosome. 
The binding sites are in data.frames with two fields, pos and weight, to indicate the binding location and strength.
}
\seealso{
\code{\link{DBChIP}}
}
\examples{
data("PHA4")
bs.list <- read.binding.site.list(binding.site.list)
str(bs.list)

}
