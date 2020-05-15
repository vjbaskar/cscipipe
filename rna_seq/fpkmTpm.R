#!/usr/bin/env Rscript
library(getopt)
mconf_installdir=Sys.getenv("mconf_installdir")
source(paste0(mconf_installdir,"/src/basic_functions.R"))

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
            Required Args:
            ----------------
            title: 
            gtfFile: gtf file using which the count table was generated
            inputFile: count table
            sampleFile: file containing <name> <condn> <library type>
            
            Input information:
            -------------------
            >> SampleFile:  <name> <condn> <library type>
            BCOR_shRNA1.star.ReadsPerGene.out.tab	BCOR_shRNA1	BCOR_shRNA	paired-end
            BCOR_shRNA2.star.ReadsPerGene.out.tab	BCOR_shRNA2	BCOR_shRNA	paired-end
            Control_shRNA1.star.ReadsPerGene.out.tab	Control_shRNA1	Control	paired-end
            Control_shRNA2.star.ReadsPerGene.out.tab	Control_shRNA2	Control	paired-end
            
            Note: <name> will not be used in this case.
            ")
    quit(status = 2)
}

title = Sys.getenv("title")
gtfFile=Sys.getenv("gtf")
countTable=Sys.getenv("inputFile")
sampleMatrix = Sys.getenv("sampleFile")

##### deseq FPKM values
message("Loading libraries")
suppressPackageStartupMessages({
    library(DESeq2)
    library(rtracklayer)
    library(GenomicFeatures)
    library("preprocessCore")
    library(ggplot2)
    library(ggthemes)
    library(reshape2)
})


# gtfFile = "~/Volumes/share/reference/gtfs/grcm38_mm10.gencode.vM16.annotation.gtf"
# countTable = "~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_RNASeq/STAR/counts.rev.frFS.Illum.tsv"
# sampleMatrix = "~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_RNASeq/STAR/samples.counttable.txt"
# 
# 
# args = commandArgs(trailingOnly = T)
# gtfFile=as.character(args[1])
# countTable = as.character(args[2])
# sampleMatrix=as.character(args[3])


" Read in GTF file "
txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
genesToConsider <- names(ebg)


" Read contrast matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
contrast.matrix = read.table(sampleMatrix)
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
x = read.table(countTable, head=1, row.names = 1)

"- Check if all the samples in the sample table are in the count table"
sdiff = setdiff(samples, colnames(x))
if(length(sdiff) > 0){
    message(" [Error] Some of the samples are not in the count table. Exiting ...")
    quit(save = "no")
    
}

counttable <- x[samples]
#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table
"--> Filter out genes that you can not obtain exon size information"
counttable <- counttable [ genesToConsider, ]
dim(counttable)
write.table(counttable, file="counts.fpkmTpm.txt", sep="\t", quote=F)

" Create DESeq object"
condition = factor(contrast.matrix$condition)
dds <- DESeqDataSetFromMatrix(countData = counttable, DataFrame(condition), design = ~ condition)
" Calculating fpkm "
" Get gene length"
" -- Including exon granges in rowRanges. You can also add in the actual gene sizes in mcols(dds)$basepairs"
temp <- rowRanges(dds)
rangesToShoveIn <- ebg[names(temp)] 
rowRanges(dds) <- rangesToShoveIn

" Normalising for size "
dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
norm_counts = counts(dds, normalized = TRUE)
deseq2_metadata = mcols(dds)


fpkmVals <- fpkm(dds)
log2fpkmVals <- log2(fpkmVals + 0.001)
qnormFpkmVals <- normalize.quantiles(as.matrix(log2fpkmVals))
rownames(qnormFpkmVals) <- rownames(log2fpkmVals)
colnames(qnormFpkmVals) <- paste0("qnorm_",colnames(log2fpkmVals))


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


fpkmMatrix <- cbind(fpkmVals, qnormFpkmVals)

" Map IDs to names"
x = import.gff2(gtfFile)
x = as.data.frame(x)
x = unique(x[,c("gene_id", "gene_name")])

" Computing tpms"
geneCounts = counts(dds)
"- getting gene lengths"
temp = rowRanges(dds)
geneWidths = sum(width(reduce(temp)))/1000

transCounts = sweep(geneCounts, 1, geneWidths, "/")
transSizeFactor = colSums(transCounts)/10^6
tpm = sweep(transCounts, 2, transSizeFactor, "/")



" Write out Raw fpkm values"
for(i in c("fpkmVals", "qnormFpkmVals", "fpkmMatrix","norm_counts", "tpm")){
    cat(i,"...\n")
    temp <- merge(x, round(get(i),3), by.x = "gene_id", by.y="row.names", all.y=T)
    write.csv(temp, file=paste0(i,".fpkm.csv"), row.names=F)
}







