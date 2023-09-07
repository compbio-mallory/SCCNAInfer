library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(optparse)
# generic script to run scope
# Define the options
option_list <- list(
  make_option(c("-b", "--bamdir"), type="character", help="bam file directory", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default="hg19", help="reference type hg38 or hg19", metavar="character"),
	make_option(c("-d", "--wd"), type="character", default="./", help="working directory", metavar="character"),
make_option(c("-o", "--out"), type="character", default="", help="output prefix", metavar="character"),
	make_option(c("-p", "--paired"), type="logical", default=FALSE, help="Set to paired-end if provided.")
)


# Parse the options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
bamPath = opt$bamdir
ref = opt$ref
paired = opt$paired
wd = opt$wd
out = opt$out

bamFile <- list.files(bamPath, pattern = '*.dedup.bam$')
bamdir <- file.path(bamPath, bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)

bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = ref)

ref_raw <- bambedObj$ref


# mapp, gc
if(file.exists("mapp.rds")){
	mapp <- readRDS("mapp.rds")
}else{
mapp <- get_mapp(ref_raw, hgref = ref)
saveRDS(mapp, "mapp.rds")
}

if(file.exists("gc.rds")){ 
  gc <- readRDS("gc.rds")
}else{
gc <- get_gc(ref_raw, hgref = ref)
saveRDS(gc, "gc.rds")
}



values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))

# converage
if(file.exists("coverageObj.rds")){
	coverageObj <- readRDS("coverageObj.rds")
}else{
	seq = "single-end"
	if(paired)
	{
		seq = "paired-end"
	}
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, 
                                  seq = seq, hgref = ref)
saveRDS(coverageObj, "coverageObj.rds")
}
Y_raw <- coverageObj$Y #read depth

# qc
if(file.exists("qcObj.rds")){
qcObj <- readRDS("qcObj.rds")
}else{
QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = Y_raw,
                    sampname_raw = sampname_raw, ref_raw = ref_raw,
                    QCmetric_raw = QCmetric_raw)
saveRDS(qcObj, "qcObj.rds")
}
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric

# get gini coefficient for each cell
Gini <- get_gini(Y)
normObj <- normalize_codex2_ns_noK(Y_qc = Y,
                                       gc_qc = ref$gc,
                                       norm_index = which(Gini<=0.12))

ploidy <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref)

normObj.scope <- normalize_scope_foreach(Y_qc = Y, gc_qc = ref$gc,
                                             K = 1, ploidyInt = ploidy,
                                             norm_index = which(Gini<=0.12), T = 1:5,
                                             beta0 = normObj$beta.hat, nCores = 8)

Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
fGC.hat <- normObj.scope$fGC.hat[[which.max(normObj.scope$BIC)]]

# cross sample
chrs <- unique(as.character(seqnames(ref)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs

for (chri in chrs) {
  message('\n', chri, '\n')
  segment_cs[[chri]] <- segment_CBScs(Y = Y,
                                      Yhat = Yhat,
                                      sampname = colnames(Y),
                                      ref = ref,
                                      chr = chri,
                                      mode = "integer", max.ns = 1)
}
iCN_sim <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
df <- as.data.frame(ref)
gc_map <- df[, c("seqnames", "start", "end", "gc", "mapp")]
gc_map$gc <- df$gc / 100
colnames(gc_map) <- c("CHROM", "START", "END", "gc", "map")
write.table(gc_map, "gc_map.tsv", sep = "\t", row.names = F, quote = F)

# select columns for chromosome, start and end
df <- df[, c("seqnames", "start", "end")]
colnames(df) <- c("CHROM", "START", "END")
CN <- cbind(df, iCN_sim)
write.table(CN, paste(wd, out, "_SCOPE_cnv.tsv", sep = "\t", row.names = FALSE, quote = F)

reads <- cbind(df, Y)
names(reads) <- names(CN)
write.table(reads, paste(wd, out, "_SCOPE_read.tsv"), sep = "\t", row.names = F, quote = F)

