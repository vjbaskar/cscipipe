#####
# Rscript deseqFPKM.R gtfFile  contrast.matrix
#####

##### deseq FPKM values
message("Loading libraries")
library(DESeq2)
library(rtracklayer)
library(GenomicFeatures)
library("preprocessCore")
library(ggplot2)
library(ggthemes)
library(reshape2)

# gtfFile = "~/Volumes/share/reference/gtfs/grcm38_mm10.gencode.vM16.annotation.gtf"
# countTable = "~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_RNASeq/STAR/counts.rev.frFS.Illum.tsv"
# contrastMatrix = "~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_RNASeq/STAR/samples.counttable.txt"


args = commandArgs(trailingOnly = T)
gtfFile=as.character(args[1])
#listOfFiles=as.character(args[2])
contrastMatrix=as.character(args[2])


" Read in GTF file "
txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
genesToConsider <- names(ebg)


" Read contrast matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
contrast.matrix = read.table(contrastMatrix)
colnames(contrast.matrix) = c("bamfiles", "names", "condition", "libtype")
rownames(contrast.matrix) = contrast.matrix$names
contrast.matrix = contrast.matrix [, c("condition","libtype") ]
head(contrast.matrix)
" conditions "
conditions = as.character(unique((contrast.matrix$condition)))
message("[ Info ] Condtions =  ", paste0(conditions, sep=" "))



" Read count table "
samples = rownames(contrast.matrix)
fileList = paste0(samples,".deseq.counts")

x = read.table(fileList[1])
rownames(x) <- x[,1]
for (i in 2:length(fileList)){
    temp = read.table(fileList[i])
    x = cbind(x, temp[,2])
}
x = x [, -1]
colnames(x) <- samples
counttable <- x
write.table(counttable, file="counts.txt", sep="\t", quote=F)

#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table
"--> Filter out genes that you can not obtain exon size information"
counttable <- counttable [ genesToConsider, ]
dim(counttable)


" -- Write out essential stats"
million <- 10^6
df <- rbind(colSums(genes.table), colSums(nonGene.table))
rownames(df) <- c("Reads to Genes", "Reads that are not mapped")
df <- melt(df)
colnames(df) <- c("Type", "Experiment", "RPM")
df$RPM <- df$RPM/million
pdf("deseq_readAlignmentStats.pdf", nrow(df)/2,4)
p <- ggplot(data = df) + aes(x = Experiment, y = RPM, group= Type ) + geom_bar(stat = "identity", aes(fill=Type), colour="black") 
p + ggtitle("Reads (per million) aligned to genes") + theme_gdocs() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + scale_fill_solarized()
dev.off()
rm(list = c("genes.table", "nonGene.table"))



" Filter count data for low count data"
countdata <- counttable
rmax = apply(countdata, 1, max)
message("[ Info ] Total genes = ", nrow(countdata))
countdata = countdata [ rmax >= 1, ]
message("[ Info ] Genes with reads = ", nrow(countdata))


" Create DESeq object"
condition = factor(contrast.matrix$condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, DataFrame(condition), design = ~ condition)
" Perform DESeq "
dds <- DESeq(dds)
dds.bkup <- dds

" Get gene length"
temp <- rowRanges(dds)
rangesToShoveIn <- ebg[names(temp)] 
rowRanges(dds) <- rangesToShoveIn

" Get FPKM "
fpkmVals <- fpkm(dds)
log2fpkmVals <- log2(fpkmVals + 0.001)
qnormFpkmVals <- normalize.quantiles(as.matrix(log2fpkmVals))
rownames(qnormFpkmVals) <- rownames(log2fpkmVals)
colnames(qnormFpkmVals) <- paste0("qnorm_",colnames(log2fpkmVals))

fpkmMatrix <- cbind(fpkmVals, qnormFpkmVals)

" Map IDs to names"
x = import.gff2(gtfFile)
x = as.data.frame(x)
x = unique(x[,c("gene_id", "gene_name")])

" Write out Raw fpkm values"
fpkm.new <- merge(x, fpkmMatrix, by.x = "gene_id", by.y="row.names", all.y=T)
dim(fpkm.new)
write.csv(fpkm.new, file="fpkm.csv", row.names=F)


" Saving Image"
save.image("deseq_fpkm.RData")


" Plot correlation "
plotcorr <- function(x){
    rn <- x$gene_id
    cn <- colnames(x)
    cn <- cn [ grepl("qnorm", cn)]
    cn
    x <- x[cn]
    rownames(x) <- rn
    temp = cor(x)
    pdf("fpkm_correlation.pdf", 5,5)
    library(gplots)
    heatmap.2(temp, trace="none", margins = c(10,20), density.info="none", main = "Pearson's correlation of expression", col=rev(redblue(20)), cexRow = 0.2 + 1/log10(10*nrow(temp)), cexCol = 0.2 + 1/log10(10*ncol(temp)))
    dev.off()
}


