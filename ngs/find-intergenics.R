library(uuutils)
library(ngsutils)
library(rtracklayer)

### create BED file of the intergenic regions (yeast-specific)

## url <- "ftp://ftp.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff"
## file <- "/data/local/genomes/yeast/cerevisiae-nofasta.gff"
file <- "~/isi.gen/public/external_data/SGD/cerevisiae_nofasta.gff"
all.gff <- import.gff3(file, genome="S.cerevisiae")

chromos <- all.gff[ all.gff$type=='chromosome',"ID"]
chr.lengths <- width(chromos)
names(chr.lengths) <- seqnames(chromos)
rm(chromos)
seqlengths(all.gff) <- chr.lengths

### get rid of all but the non-dubious protein_codings
subset <- (all.gff$type == "gene" & all.gff$orf_classification != "Dubious")
gff <- all.gff[subset,]
names(gff) <- values(gff)$ID
gff[ is.na(gff$gene) ]$gene <- gff[ is.na(gff$gene) ]$ID # just systematic there

### --- create bed file with intergenic regions, and name
### --- them by flanking genes and type (convergent, divergent, tandem)
genes.nostrand <- gff
strand(genes.nostrand) <- '*'
interg <- gaps(genes.nostrand)
interg <- interg[ strand(interg)=="*" ] # also returns useless "+" and "-" 

values(interg) <- DataFrame(left="", left.strand="",
                            right="", right.strand="",
                            name="", type="")

f <- follow(interg, genes.nostrand)
interg[!is.na(f)]$left = gff[ f[!is.na(f)] ]$gene
interg[!is.na(f)]$left.strand = strand(gff[ f[!is.na(f)] ])

p <- precede(interg, genes.nostrand)
interg[!is.na(p)]$right = gff[ p[!is.na(p)] ]$gene
interg[!is.na(p)]$right.strand = strand(gff[ p[!is.na(p)] ])

interg[ interg$left.strand=="" | interg$right.strand==""  ]$type <-
  'telomeric'

interg[ interg$left.strand == interg$right.strand ]$type <-
  'tandem'

interg[ interg$left.strand == '-' &
       interg$right.strand == '+'  ]$type <-  'divergent'

interg[ interg$left.strand == '+' &
       interg$right.strand == '-'  ]$type <-  'convergent'

makename <- function(left,right,sep)paste0(left,sep,right)

interg[ interg$type=="telomeric"  ]$name <-
  with(values(interg[ interg$type=="telomeric" ]), makename(left, right, sep="*"))

interg[ interg$type=="divergent"  ]$name <-
  with(values(interg[ interg$type=="divergent" ]), makename(left, right, sep="<>"))

interg[ interg$type=="convergent"  ]$name <-
  with(values(interg[ interg$type=="convergent" ]), makename(left, right, sep="><"))

subset <- interg$type=="tandem" &  interg$left.strand=="+"
interg[ subset ]$name <-
  with(values(interg[ subset ]), makename(left, right, sep=">>"))
subset <- interg$type=="tandem" &  interg$left.strand=="-"
interg[ subset ]$name <-
  with(values(interg[ subset ]), makename(left, right, sep="<<"))

export.bed(con="intergenic2.bed", interg)

### The next bit really should also be rewritten as for the
### intergenic thing! 100x faster and cleaner ...

## allocate storage for properties
all.genes <- names(gff)
n <- length(all.genes)

### 5'-related info:
vP.name <- rep("", n); names(vP.name) <- all.genes
vP.type <- rep("", n); names(vP.type) <- all.genes
vP.dist <- rep(-1, n); names(vP.dist) <- all.genes
vP.overlap <- rep("", n); names(vP.overlap) <- all.genes

### 3'-related info:
iiiP.name <- rep("", n); names(iiiP.name) <- all.genes
iiiP.type <- rep("", n); names(iiiP.type) <- all.genes
iiiP.dist <- rep(-1, n); names(iiiP.dist) <- all.genes
iiiP.overlap <- rep("", n); names(iiiP.overlap) <- all.genes
rm(n)

chr.names <- as.vector(sort(unique(seqnames(all.gff))))
for (chr.name in chr.names) {
    chr <- gff[seqnames(gff)==chr.name]
    gene.names <- names(chr)

    for(gene.name in gene.names) {
        gene <- chr[gene.name]
        strand <- as.vector(strand(gene))
        stopifnot(strand %in% c('+', '-'))
        fwd <- strand=='+'

        ## check neighbours at 5' end:
        vP <- resize(gene, width=1, fix='start') # start means 5' here !!@#$
        names(vP) <- gene.name
        idx <- (if(fwd)follow else precede)(vP, chr, ignore.strand=TRUE)
        if(!is.na(idx)) {                   #otherwise: end of chromosome
            vP.ngb <- chr[idx]                #closest neighbouring gene 
            vP.name[gene.name] <- names(vP.ngb)
            vP.dist[gene.name] <- distance(vP, vP.ngb, ignore.strand=TRUE)
            vP.overlap[gene.name] <- (length(findOverlaps(vP, chr, ignore.strand=TRUE)) >= 2)
            vP.type[gene.name] <- ifelse(strand(vP.ngb)==strand(vP), "3",  "5")
### type means: what kind of "end" (5' or 3') the gene closest to
### the current gene's 5'-end has. E.g if that is 5', we have a divergent 
### promoter.
        }

        ## same for neighbours at 3'-end
        iiiP <- resize(gene, width=1, fix='end') # start means 3' here !!@#$
        names(iiiP) <- gene.name
        idx <- follow(iiiP, chr)
        idx <- (if(fwd)precede else follow)(iiiP, chr, ignore.strand=TRUE)
        if(! is.na(idx)) {
            iiiP.ngb <- chr[idx]
            iiiP.name[gene.name] <- names(iiiP.ngb)
            iiiP.dist[gene.name] <- distance(iiiP, iiiP.ngb, ignore.strand=TRUE)
            iiiP.overlap[gene.name]<- (length(findOverlaps(iiiP, chr, ignore.strand=TRUE))>=2)
            iiiP.type[gene.name] <- ifelse(strand(iiiP.ngb)==strand(iiiP), "5",  "3")
            ## as before, but now "end" of the gene closest to the current gene's
            ## 3'-end. E.g. if it's 3', we have a convergent intergenic region.
        }
    }
}                                       #for chr.name

vP.dist[is.na(vP.dist)] <- -1 
iiiP.dist[is.na(iiiP.dist)] <- -1 

geneSymbol <- gff$gene
missing <- is.na(geneSymbol)
geneSymbol[missing] <- names(gff)[missing]
names(geneSymbol) <- names(gff)

df <- data.frame(systematicName=all.genes,
                 geneSymbol=geneSymbol,
                 chr=seqnames(gff),
                 strand=strand(gff),
                 length=width(gff),
                 vP.dist=vP.dist,
                 vP.type=vP.type,
                 vP.overlap=vP.overlap,
                 vP.geneSymbol=geneSymbol[vP.name],
                 iiiP.dist=iiiP.dist,
                 iiiP.type=iiiP.type,
                 iiiP.overlap=iiiP.overlap,
                 iiiP.geneSymbol=geneSymbol[iiiP.name],
                 desc=values(gff)$display,
                 row.names=all.genes)

## write.tab(file="/data/local/genomes/yeast/promoter-data.txt", df)
write.tab(file="promoter-data.txt", df)


