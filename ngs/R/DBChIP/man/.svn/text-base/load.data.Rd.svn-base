\name{load.data}
\alias{load.data}
\title{
load.data
}
\description{
Assemble data for further analysis.
}
\usage{
load.data(chip.data.list, conds, consensus.site, input.data.list = NULL, data.type = "MCS", chr.vec = NULL, chr.exclusion = NULL, chr.len.vec = NULL, norm.factor.vec = NULL, frag.len = 200)
}
\arguments{
  \item{chip.data.list}{
a list of ChIP data where each list item corresponds to one ChIP library. The name of the items should be unique. Biological replicates should be in separate items. Each item can be one of three accepted data types: MCS, AlignedRead and BED.
}
  \item{conds}{
a vector of conditions of ChIP libraries. Should be the same order as chip.data.list, or the names should be specified as a permutation of the names of chip.data.list.
}
  \item{consensus.site}{
consensus binding sites. Should be the result of \link{site.merge}.
}
  \item{input.data.list}{
a list of control data. Should have same data type as in chip.data.list. The names of the items should be unique, and each name should be matched to either a ChIP replicate name when the two are paired or a condition name in general.
}
  \item{data.type}{
"MCS", "AlignedRead" or "BED". See Details.
}
  \item{chr.vec}{
a vector of chromosomes in data. User can specify chr.vec, or it can be computed from the ChIP and control samples.
}
  \item{chr.exclusion}{
user can either specify chr.vec, or specify the chromosomes to exclude through this parameter.
}
  \item{chr.len.vec}{
a vector of chromosome lengths corresponding to chr.vec. Can be specified if known, or will be computed as the largest 5' end position in the data.
}
  \item{norm.factor.vec}{
a vector of normalization factors between the ChIP and control libraries when controls are available. Can be specified, or will be computed by DBChIP.
}
  \item{frag.len}{
average fragment length. Default 200 bp.
}
}
\details{
The ChIP and control data should be properly filtered before the analysis to avoid artifacts.
For example, reads mapping to mitochondrial DNA, or Y chromosome for female samples will need to be filtered.\cr

Filtering of chromosomes can be achieved through specification of \code{chr.vec} and/or \code{chr.exclusion}. Only reads from chromosomes in \code{chr.vec} but not in \code{chr.exclusion} are utilized in the analysis.\cr

User can include or exclude sex chromosomes in the computation, depending on whether protein-DNA bindings on sex chromosomes are of research interest.\cr

Biological replicates of a ChIP sample should be kept separate so that dispersion can be properly estimated. On the other hand, replicates of a control/input sample should be merged because the purpose of the control samples is to estimate the background for testing and plotting. 
One exception would be when a control replicate is paired with a ChIP replicate, for example, they are coming from the same batch, a portion of which is used for IP and the other portion is used for control. In such case, the control replicate can be kept separate with the same name of the matching ChIP replicate.

data.type
\itemize{
\item \code{MCS} Minimum ChIP-Seq format. data.frame with fields: chr (factor), pos (integer) and strand (factor, "+" and "-"). pos is 5' location. This is different from eland default which use 3' location for reverse strand.
\item \code{AlignedRead} from Bioconductor ShortRead package (with support of commonly used formats, including Eland, MAQ, Bowtie, SOAP and BAM).
\item \code{BED} with at least first 6 fields (chrom, start, end, name, score and strand), http://genome.ucsc.edu/FAQ/FAQformat.html#format1.
}
}
\value{
A list with following components:
\item{chip.list }{the list of ChIP data in internal MCS format.}
\item{conds }{vector of conditions of ChIP libraries.}
  \item{frag.len}{
average fragment length.
}
  \item{chr.vec}{
a vector of chromosomes in data.
}
  \item{chr.len.vec}{
a vector of chromosome lengths corresponding to chr.vec.
}
  \item{consensus.site}{
consensus sites. It is a data.frame, for details, see \code{site.merge}.
}
  \item{input.list}{
the list of control data. The components from this and below are only available when control data are available.
}
  \item{matching.input.names}{
the matching input names for ChIP replicates.
}
  \item{chip.background.size}{
the background sizes (number of reads) for ChIP replicates. Background is defined by exluding the neighborhoods of the consensus sites.
}
  \item{input.background.size}{
the background sizes for input replicates.
}
  \item{norm.factor.vec}{
vector of normalization factors between the ChIP and control libraries when controls are available. Can be specified, or will be computed by DBChIP.
}

}

\seealso{
\code{\link{DBChIP}}.
}

\examples{
data("PHA4")
conds <- factor(c("emb","emb","L1", "L1"), levels=c("emb", "L1"))

bs.list <- read.binding.site.list(binding.site.list)

## compute consensus site
consensus.site <- site.merge(bs.list, in.distance=100, out.distance=250)

#load data
dat <- load.data(chip.data.list=chip.data.list, conds=conds, consensus.site=consensus.site, input.data.list=input.data.list, data.type="MCS")
names(dat)
}

\keyword{ load.data }

