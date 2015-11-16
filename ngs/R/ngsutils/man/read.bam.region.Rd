\name{read.bam.region}
\alias{read.bam.region}
\title{
  Read a subset of a bam file limited to a particular genomic region.
}

\description{
  Reading a complete bam file often costs too much memory. This function
  you can read a subset, typically a chromosome or chromosome arm. The
  bam file must be indexed, i.e., there must be a .bam.bai file with the
  same name.
}

\usage{
  read.bam.region(file, region = NULL)
}

\arguments{
  \item{file}{
    The name of the bam file. The bam file must be index,
  i.e. file.bam.bai must also exist.
}
  \item{region}{

    One or more genome regions, specified as for samtools. Commas to
  separate thousands are allowed. If region is \code{NULL}, \code{NA} or
  the empty string, the whole file is read. If region contains just the
  chromosome name, the whole chromosome is read. If \code{region} is of
  length > 1, they must all have the same 'shape'. E.g., c("chrX",
  "chrXI:100,000-200,000") is invalid.  } } \value{ A
  \code{\link{GRanges}} object is returned. Only seqname, start, end and
  strand are returned, there is no metadata.
}

\author{
 Philip Lijnzaad <plijnzaad@gmail.com>
}

\examples{
\dontrun{
  gr <- read.bam.region(file="foo.bam", region=c("chrIX:10,000-20,000",
                                                 "chrX:10,000-20,000"))
}
}

\keyword{misc}

