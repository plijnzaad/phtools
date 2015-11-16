\name{DBChIP}
\alias{DBChIP}
\title{ DBChIP }
\description{
Detecting differential binding of transcription factors with ChIP-seq
}
\usage{
DBChIP(binding.site.list, chip.data.list, conds, input.data.list = NULL,
 data.type = c("MCS", "AlignedRead", "BED"), frag.len = 200, chr.vec = NULL, 
 chr.exclusion = NULL, chr.len.vec = NULL, subtract.input = FALSE, norm.factor.vec = NULL, 
 in.distance = 100, out.distance = 250, window.size = 250,
 dispersion=NULL, common.disp=TRUE, prior.n=10, 
 two.sample.method="composite.null", allowable.FC=1.5, collapsed.quant=0.5)
}
\arguments{
  \item{binding.site.list}{
a list of data.frames. Each data.frame corresponds to one condition in comparison 
and has three fields, chr, pos, and weight, to indicate the binding location and strength.
}
  \item{chip.data.list}{
a list of ChIP data where each list item corresponds to one ChIP library. The name of the items should be unique. Biological replicates should be in separate items. Each item can be one of three accepted data types: MCS, AlignedRead and BED.
}
  \item{conds}{
a vector of conditions of ChIP libraries. Should be the same order as chip.data.list, or the names should be specified as a permutation of the names of chip.data.list.
}
  \item{input.data.list}{
a list of control data. Should have same data type as in chip.data.list. The names of the items should be unique, and each name should be matched to either a ChIP replicate name when the two are paired or a condition name in general.
}
  \item{data.type}{
  "MCS", "AlignedRead" or "BED". See Details.
}
  \item{frag.len}{
average fragment length. Default 200 bp.
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
  \item{subtract.input}{
logical. Whether input will be subtracted from ChIP when counting the binding reads. Default is FALSE.
}
  \item{norm.factor.vec}{
a vector of normalization factors between the ChIP and control libraries when controls are available. Can be specified, or will be computed by DBChIP.
}
  \item{in.distance}{
within cluster distance. If the distance between centroids of two clusters are smaller than this value, the clusters will be merged into one. Default value 100 bp.
}
  \item{out.distance}{
out of cluster distance. If the distance between centroids of two clusters are larger than this value, they are considered different binding sites. Also double-function as the window size to count reads around each site. Default value 250 bp.
}
  \item{window.size}{
the window size to count reads around a binding site.
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

Users are recommended to study the histogram of the $p$-values for model checking. More specifically, the $p$-values between 0.5 and 1 should be roughly uniform. 
When many replicates are available, users can also randomly split biological replicates of the same condition and perform comparisons through DBChIP using the estimated dispersion parameter to check whether the $p$-values look uniform.

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
  \item{site.count}{
a data.frame of read counts at each consensus site.
}
  \item{test.stat}{
a data.frame of test statistics for testing non-differential binding at each site, include p-values and fold changes.
}
  \item{input.list}{
the list of control data. The components from this and below are only available when control data are available.
}
  \item{matching.input.names}{
the matching input names for ChIP libraries.
}
  \item{norm.factor.vec}{
vector of normalization factors between the ChIP and control libraries when controls are available. Can be specified, or will be computed by DBChIP.
}

}
\references{
Liang, K and Keles, S (2012). Detecting differential binding of transcription factors with ChIP-seq. \emph{Bioinformatics}, 28, 121-122.
}
\author{
Kun Liang, \email{kliang@stat.wisc.edu}
}

\seealso{
\code{\link{DBChIP-package}}, \code{\link{PHA4}}, \code{\link{read.binding.site.list}}, \code{\link{site.merge}}.
}

\examples{
data("PHA4")
dat <- DBChIP(binding.site.list, chip.data.list=chip.data.list, input.data.list=input.data.list, conds=conds, data.type="MCS")
rept <- report.peak(dat)
rept
#pdf("Diff.Binding.pdf")
plotPeak(rept, dat)
#dev.off()

## experienced users can proceed in a step by step fashion such that if program
## needs to be run for a different setting, intermediate results can be saved and reused.
data("PHA4")
conds <- factor(c("emb","emb","L1", "L1"), levels=c("emb", "L1"))

bs.list <- read.binding.site.list(binding.site.list)

## compute consensus site
consensus.site <- site.merge(bs.list, in.distance=100, out.distance=250)

dat <- load.data(chip.data.list=chip.data.list, conds=conds, consensus.site=consensus.site, input.data.list=input.data.list, data.type="MCS")

## count ChIP reads around each binding site
dat <- get.site.count(dat, window.size=250)

## test for differential binding
dat <- test.diff.binding(dat)

# report test result and plot the coverage profiles
rept <- report.peak(dat)
rept
plotPeak(rept, dat)
}
