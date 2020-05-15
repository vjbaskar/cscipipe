
#####
# Rscript deseq.R gtfFile countTable.txt contrast.matrix condition1 condition2
#####


#!/usr/bin/env Rscript



library(getopt)

args_return=Sys.getenv("args_return")
if(args_return > 0){
		message("Help: 
		title: 
		gtfFile: gtf file using which the count table was generated
		inputFile: count table
		sampleFile: file containing <name> <condn> <library type>
		
		")
		quit(status = 2)
}

title = Sys.getenv("title")
message(title)
gtfFile=Sys.getenv("gtf")
countTable=Sys.getenv("inputFile")
contrastMatrix = Sys.getenv("sampleFile")




##### deseq FPKM values
message("Loading libraries")
suppressMessages({
	library(DESeq2)
	library(rtracklayer)
	library(GenomicFeatures)
	library("preprocessCore")
	library(ggplot2)
	library(ggthemes)
	library(reshape2)
})

#gtfFile = "~/nfs_vijay/GENOMES/SANGER/GRCm38/ensembl_84_transcriptome-GRCm38.gtf"
#countTable = "counts.txt"
#contrastMatrix = "samples.GB_Imran.deseq.txt"
#c1="brainControl"
#c2="brainTumour"


" Read in GTF file "
gtfData <- import.gff2(gtfFile)
annot_data <- as.data.frame(gtfData)
annot_data <- annot_data[ annot_data$type == "gene", ]
head(annot_data)
annot_data <- data.frame("gene_name" = annot_data$gene_name,
                        "gene_id" = annot_data$gene_id,
                        "coords" = paste0(annot_data$seqnames, ":", annot_data$start, "-", annot_data$end), 
                        "strand" = annot_data$strand,
                        "source" = annot_data$source
                        )
annot_data = unique(annot_data)
genesToConsider = as.character(annot_data$gene_id)

" Read count table "
x = read.table(countTable, head=T, row.names=1)
#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table
"--> Filter out genes that you can not obtain exon size information"
counttable <- counttable [ genesToConsider, ]
dim(counttable)

" Read contrast matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
contrast.matrix = read.table(contrastMatrix)
colnames(contrast.matrix) = c("names", "condition", "libtype")
rownames(contrast.matrix) = contrast.matrix$names
contrast.matrix = contrast.matrix [, c("condition","libtype") ]
#head(contrast.matrix)
" conditions "
conditions = as.character(unique((contrast.matrix$condition)))
message("[ Info ] Conditions =  ", paste0(conditions, sep=" "))


" Extract data corresponding to the conditions"
colData = subset(contrast.matrix, condition %in% c(c1,c2))
head(colData)

" Filter count data for low count data, reads/gene > 0"
countdata <- counttable
countdata = as.matrix(countdata [ , rownames(colData)])

rmax = apply(countdata, 1, max)
message("[ Info ] Total genes = ", nrow(countdata))
countdata = countdata [ rmax >= 1, ]
message("[ Info ] Genes with reads = ", nrow(countdata))


" -- Check if the row order in contrast and columns in count data are in the same order"
if( ! all(rownames(colData) == colnames(countdata)) ){
    message("[ Error ] Total rows in filtered contrast matrix != Total columns in count table")
    message("[ Error ] Aborting")
    save.image("deseq_error.RData")
    quite(save="no")
} else {
    message("[ Info ] All seems correct. Proceeding forward ...")
}

" Create DESeq object"
condition = factor(colData$condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~ condition)
dds$condition <- factor(dds$condition, levels=c(c1,c2))
" Perform DESeq "
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
res.bkup <- res
#summary(res)




" Extract results"
genesInRes = as.character(rownames(res))
" -- Read FPKM data"
fpkmMatrix <- read.csv(file = "fpkm.csv")
fpkmMatrix <- fpkmMatrix[ fpkmMatrix$gene_id %in% genesInRes,  ]
fpkmMatrix <- fpkmMatrix [ , ! grepl("qnorm", colnames(fpkmMatrix)) ]
dim(fpkmMatrix)
fpkmMatrix_metadata = fpkmMatrix [,c(1:2)]
cn1 <- rownames(colData [ colData$condition == c1, ])
cn2 <- rownames(colData [ colData$condition == c2, ])

med.FPKM.c1 <- rowMedians(as.matrix(fpkmMatrix[,cn1]))
med.FPKM.c2 <- rowMedians(as.matrix(fpkmMatrix[,cn2]))

medFPKMmatrix = cbind(fpkmMatrix_metadata, cbind(med.FPKM.c1, med.FPKM.c2))
colnames(medFPKMmatrix) <- c(colnames(fpkmMatrix_metadata), c1, c2)


res.filt <- function(res.data, fc, p.adj, annot_data, fpkmMatrix) {
    temp = as.data.frame(subset(res.data, abs(log2FoldChange) >= fc & padj <= p.adj ))
    x = rownames(temp)
    x = cbind(x, temp)
    colnames(x) = c("gene", colnames(temp))
    x = as.data.frame(x)
    x = merge(annot_data, x, by.x = "gene_id",by.y="gene", all.y = T)
    fpkmMatrix <- fpkmMatrix[, ! grepl("gene_name", colnames(fpkmMatrix)) ]
    x = merge(x, fpkmMatrix, by="gene_id", all.x=T)
    x = x [ order(x$log2FoldChange, x$padj), ]
    return(x)
}

res.lfcThreshold <- function(dds, fc, p.adj, annot_data, fpkmMatrix) {
    res.data <- results(dds, lfcThreshold = fc, alpha = p.adj, format = "DataFrame")
    temp = as.data.frame(res.data)
    temp <- subset(temp, (abs(log2FoldChange) >= fc & padj <= p.adj))
    x = rownames(temp)
    x = cbind(x, temp)
    colnames(x) = c("gene", colnames(temp))
    x = as.data.frame(x)
    x = merge(annot_data, x, by.x = "gene_id",by.y="gene", all.y = T)
    fpkmMatrix <- fpkmMatrix[, ! grepl("gene_name", colnames(fpkmMatrix)) ]
    x = merge(x, fpkmMatrix, by="gene_id", all.x=T)
    x = x [ order(x$log2FoldChange, x$padj), ]
    return(x)
}



" -- Data from DESeq "
res.results = list()
res.results[["all"]] = res.filt(res, 0,1, annot_data, medFPKMmatrix)
res.results[["abs.logFC >= 1 + padj <= 0.05"]] = res.filt(res, 1, 0.05, annot_data, fpkmMatrix)
res.results[["abs.logFC >= 1 + padj <= 0.01"]]  = res.filt(res, 1, 0.01, annot_data, fpkmMatrix)
res.results[["abs.logFC >= 2 + padj <= 0.05"]]  = res.filt(res, 2, 0.05, annot_data, fpkmMatrix)
res.results[["abs.logFC >= 2 + padj <= 0.01"]]  = res.filt(res, 2, 0.01, annot_data, fpkmMatrix)
res.results[["lfcThreshold abs.logFC >= 0.5"]] = res.filt(res, 0.5, 1, annot_data, fpkmMatrix)
res.results[["lfcThreshold abs.logFC >= 1"]] = res.lfcThreshold(dds, 1, 0.99, annot_data, fpkmMatrix)
res.results[["lfcThreshold abs.logFC >= 1.5"]] = res.lfcThreshold(dds, 1.5, 0.99, annot_data, fpkmMatrix)
res.results[["lfcThreshold abs.logFC >= 2"]] = res.lfcThreshold(dds, 2, 0.99, annot_data, fpkmMatrix)
res.results[["lfcThreshold abs.logFC >= 2.5"]] = res.lfcThreshold(dds, 2.5, 0.99, annot_data, fpkmMatrix)
res.results[["lfcThreshold abs.logFC >= 3"]] = res.lfcThreshold(dds, 3, 0.99, annot_data, fpkmMatrix)




" Write results "
fold = paste(c1, "VS", c2, sep="_")
system(paste("mkdir ", fold), intern = F)
library(WriteXLS)
WriteXLS(res.results, ExcelFileName=paste(fold,"/deseqAll.xlsx",sep=""), row.names=F)

" Plots "
plot.volcano <- function(res, lfc, pval){
	tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
	head(tab)
	par(mar = c(5, 4, 4, 4))
	plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
	signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
	points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
	abline(h = -log10(pval), col = "green3", lty = 2) 
	abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
	mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
	mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
}

pdf(paste(fold,"/deseq.pdf",sep=""),5,5)
plotMA(res, frame=F, alpha=0.05)
plot.volcano(res, 1, 0.05)
plotDispEsts(dds)
dev.off()


" Save R image"
save.image(paste(fold,"/deseq.RData",sep=""))
