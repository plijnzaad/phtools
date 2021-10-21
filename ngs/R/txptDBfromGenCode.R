### Adapted from https://chitka-kalyan.blogspot.nl/2014/02/creating-gencode-transcript-database-in.html
stop("did not work:   ''Required tables tx, tx2exon, chromosome are not present in the database! ''")

library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)

###Get the chromosome info as a dataframe
###One can use the script fethChromsomeSize from UCSC to the get info, filter it to remove information from 
###non-std chromosomes
###Add the header "chrom length  is_circular" additional column to say that none of the chromosomes are circular.
###  chrom length  is_circular
###  chr1 249250621       FALSE
###  chr2 243199373       FALSE
###  chr3 198022430       FALSE
file <- "/hpc/local/CentOS7/gen/data/genomes/human/ensembl/GRCh38/dna/hg38p5.chromSizes"
chrom.info <- read.table(file=file, header=TRUE)
###Create the transcriptdb from the gencode gtf file
###Download the latest gencode comprehensive gtf file from gencode website

gtf <- '/hpc/local/CentOS7/gen/data/genomes/human/gencode/gencode.v26.annotation-sorted.gtf'
source <- 'adapted from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz'

save.as <- '/hpc/local/CentOS7/gen/data/genomes/human/gencode/hg38-gencode-v26.sqlite'

gencode<-makeTxDbFromGFF(gtf,
                         format="auto",
                         dataSource=source,
                         chrominfo=chrom.info, 
                         organism="Homo sapiens")

###Save the transcriptDb object as a sql database object
saveDb(gencode, file=save.as)
stop("till here")
loadDb(save.as)

#### Create GRAnges objects for the each feature
genes.gencode<-sort(genes(gencode))
intergenic.gencode<-gaps(genes.gencode)
transcript.gencode<-sort(transcripts(gencode))
cds.gencode<-sort(cds(gencode))
exons.gencode<-sort(exons(gencode))
introns.gencode<-sort(unlist(intronsByTranscript(gencode)))
###save the combined obj for easy loading in future.
save(genes.gencode,intergenic.gencode,transcript.gencode,cds.gencode,exons.gencode,introns.gencode, file="/home/kalyankpy/reg_elements/gencode_all_features.rda")
###Load the database in future
load("/home/kalyankpy/reg_elements/gencode_all_features.rda")
###Alternatively save all the above objects into a single 'GenomicRangesList'
annotation.all.features=GenomicRangesList(genes.gencode,intergenic.gencode,transcript.gencode,cds.gencode,exons.gencode,introns.gencode)
###save this list of features
save(annotation.all.features, file="/home/kalyankpy/reg_elements/gencode_all_features_list.rda")
###Load in future
load("/home/kalyankpy/reg_elements/gencode_all_features_list.rda")



txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                            dataset="hsapiens_gene_ensembl",
                            transcript_ids=NULL,
                            circ_seqs=DEFAULT_CIRC_SEQS,
                            filter=NULL,
                            id_prefix="ensembl_",
                            host="www.ensembl.org",
                            port=80,
                            taxonomyId=NA,
                            miRBaseBuild=NA)
     


getChromInfoFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
