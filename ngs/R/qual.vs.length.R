#!/usr/bin/env Rscript

### simple script to plot the quality vs. fragment length, including marginal histograms

library(uuutils)

file <- Sys.getenv("file")
if(file=="")
  stop("Usage: env file=filename.txt qual.vs.length.R\nFile is assumed to be tab-delimited, 1st column the MAPQ, 2nd the TLEN\n")


warning("reading ", file , " ...") 
x <- read.table(file)
names(x) <- c("qual", "length")

sample <- strip.extension(file)

f <- paste0(sample, ".pdf")
warning("creating ", f) 
pdf.a4landscape(f)

## sub-sample, otherwise too much to plot
s <- sample(nrow(x), 100000)
with(x[s,], plot.cor(x=jitter(abs(length)+10), y=jitter(qual),
                     pch=1, cex=0.7, show.p.values=FALSE, show.regression=FALSE, show.line=FALSE, confidence.interval=NA,
                     sub=paste0(sample," (subsample)"), xlab="fragment length", ylab="mapping quality"
                 # , xlim=c(0,100)
                 ))

### add stats (originally written for Bowtie2 output, may need tweaking)
n <- nrow(x)
qhigh <- with(x, sum(qual >= 40))
qmid <- with(x, sum(qual>1 & qual < 40))
q1 <- with(x, sum(qual==1))
q0 <- with(x, sum(qual==0))


p <- function(s, n,tot) {
    maxl <- nchar(commafy(tot))
    sprintf("%s: %*s ( %4.2g%% )", s, maxl, commafy(n), 100*n/tot)
}

txt <- paste(sep="\n",
               p(" q>=40", qhigh,n),
               p("2<q<40", qmid,n),
               p("  q==1", q1,n),
               p("  q==0", q0,n)
               )
cat("qualities:\n",txt)

par(family="mono")
text(x=1200,y=30, txt)
dev.off()
sink(file=stderr())
sessionInfo()
