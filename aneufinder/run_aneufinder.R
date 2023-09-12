library(AneuFinder)
    library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)
# generic script to run scope
# Define the options
option_list <- list(
  make_option(c("-b", "--bamdir"), type="character", help="bam file directory", metavar="character"),
        make_option(c("-d", "--wd"), type="character", default="./", help="working directory", metavar="character"),
  make_option(c("-w", "--binwidth"), type = "character", help = "binning variable")
)


# Parse the options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
bamPath = opt$bamdir
wd = opt$wd

var.width.ref <- system.file("extdata", "hg19_diploid.bam.bed.gz", package="AneuFinderData")
blacklist <- system.file("extdata", "blacklist-hg19.bed.gz", package="AneuFinderData")
datafolder<-bamPath
outputfolder<-wd
    Aneufinder(inputfolder = datafolder, outputfolder = outputfolder, assembly = 'hg19', numCPU = 4, binsizes = c(as.numeric(opt$binwidth)), chromosomes = c("chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"), blacklist = blacklist, correction.method = 'GC', GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19, refine.breakpoints=FALSE, method = 'edivisive')
