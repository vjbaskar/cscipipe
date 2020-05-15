#!/usr/bin/env Rscript
library(getopt)
mconf_installdir=Sys.getenv("mconf_installdir")
source(paste0(mconf_installdir,"/src/basic_functions.R"))

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
* Description:
--------------
Creates fpkm and tpm from count data.
RPKM:
----
raw counts = counts + pseudocount (default = 1)

1. Find genes with > pseudocount values. 
2. Get the q(75) of it. q75 <- quantile(vec/sum(vec), 0.75)
3. rpkm[gene,sample] = (count[gene,sample]+1)/(q75[sample]*10^6*geneLength[gene])

TPM:
---
raw counts = counts + pseudocount (default = 1)
1. length_normalised_counts = raw_counts/genelengths_kb
2. Remove genes with > 0 values.
3. q75 <- quantile(vec/sum(vec), 0.75)
4. tpm[gene,sample] = (length_normalised_counts[gene,sample]+1)/q75[sample]


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

#### For laptop
#setwd("/Volumes/scratch119/ETIENNE/2018.02.PPM1D.EtienneJoanne.RNASeq/TPMS/")
#title = "temp"
#gtfFile = "ensembl_75_transcriptome-1000Genomes_hs37d5.gtf"
#countTable = "featurecount_table.txt"
#sampleMatrix = "samples.txt"
####

message("Loading libraries")
suppressMessages({
    library(DESeq2)
    library(rtracklayer)
    library(GenomicFeatures)
    library(WriteXLS)
    library(preprocessCore)
})
# title = "trial1"
# gtfFile = "ensembl_75_transcriptome-1000Genomes_hs37d5.gtf"
# countTable = "featurecount_table.txt"
# sampleMatrix = "samples.txt"



" Read sample matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
sample.matrix = read.table(sampleMatrix)
colnames(sample.matrix) = c("bam","names", "condition", "libtype")
rownames(sample.matrix) = sample.matrix$names
sample.matrix = sample.matrix [, c("condition","libtype") ]
#head(sample.matrix)
" conditions "
conditions = as.character(unique((sample.matrix$condition)))
message("[ Info ] Conditions =  ", paste0(conditions, sep=" "))
#sample.matrix

" Read count table "
info("Reading count table")
x = read.table(countTable, head=T, row.names=1)
info("Check if your count file is read correctly")
x[1:5,1:5]
#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table [, rownames(sample.matrix) ]

" Read in GTF file "
info("Get gene lengths for normalisation .... This takes time")
gtfData <- import.gff2(gtfFile)
annot_data <- as.data.frame(gtfData)
annot_data <- annot_data[ annot_data$type == "gene", ]
annot_data <- data.frame("gene_name" = annot_data$gene_name,
                         "gene_id" = annot_data$gene_id,
                         "coords" = paste0(annot_data$seqnames, ":", annot_data$start, "-", annot_data$end), 
                         "strand" = annot_data$strand,
                         "source" = annot_data$source
)
annot_data = unique(annot_data)

" -- Convert to txdb"
txdb <- makeTxDbFromGRanges(gtfData)
" -- Convert to genes GRList"
genes_txdb <- exonsBy(txdb, by = "gene")
" -- Merge exonic regions per gene "
genes_fusetxdb <- sapply(genes_txdb, reduce)
" -- Computing length for the genes" 
geneLengths <- sapply(genes_fusetxdb, function(x) sum(width(x)))
" -- gene lengths available for: "
genesToConsider = as.character(names(geneLengths))

"--> Filter out genes that you can not obtain exon size information"
cat(" Total genes for which length available: ", length(genesToConsider),"\n")
cat(" Total genes in count table : ", nrow(counttable),"\n")
commonGenes <- genesToConsider [ genesToConsider %in% rownames(counttable)]
cat(" Total genes in common : ", length(commonGenes),"\n")
counttable <- counttable [ commonGenes, ]
geneLengths <- geneLengths [ commonGenes ]


colData <- data.frame(
    "condition" = colnames(counttable),
    "libtype" = "paired-end"
)
rownames(colData) <- colnames(counttable)
dds <- DESeqDataSetFromMatrix(countData = counttable, colData = colData, design = ~ 1)
dds <- estimateSizeFactors(dds)
geneCounts <- counts(dds, normalized = FALSE)
normCounts <- counts(dds, normalized = TRUE)


mergeWithAnnotation <- function(matx, annot_data){
    temp = rownames(matx)
    temp <- sapply(temp, function(x) strsplit(x, "[.]")[[1]][1] ) 
    rownames(matx) <- temp
    annot_data$gene_id <- sapply(as.character(annot_data$gene_id), function(x) strsplit(x, "[.]")[[1]][1] ) 
    annot_matx <- merge(annot_data, matx, by.x = "gene_id", by.y = "row.names", all.y=T)
    return(annot_matx)
}



info("Calculate TPMS")

# function for upper quartile normalisation
UQ <- function(y, cutoff){ # x is number vector
    temp <- y [ y > cutoff ]
    q75 <- quantile(temp/sum(temp), 0.75)
    return(q75)
}

pseudoCount = 1
rawCounts = geneCounts + pseudoCount
geneLengths.kb <- geneLengths [ rownames(rawCounts) ] / 10^3 # geneLengths in kb
tpm <- rawCounts/geneLengths.kb  # Normalise for gene lengths
sf <- apply(tpm, 2, cutoff = 0, UQ)

dataNorm <- sf # The column/expt wise normalised is on the geneLength-normalised-matrix

#tpm <- tpm/dataNorm
for(i in colnames(tpm)){
    message("Normalising for ", i)
    tpm[,i] = tpm[,i]/(dataNorm[i]*10^6)
}
tpm_annot = mergeWithAnnotation(tpm, annot_data) # merge tpm=rownames and annot_data = gene_id
info(" Write TPM to csv")
write.table(tpm, file=paste(title,".tpm.tsv",sep=""), sep="\t", quote=F)
write.csv(tpm, file=paste(title,".tpm.csv",sep=""))
write.csv(tpm_annot, file=paste(title,"annot.tpm.csv",sep=""), row.names=F)

info("Calculating FPKMs")
dataNorm <- apply(rawCounts, 2,cutoff = pseudoCount, UQ)
rpkm <- rawCounts/geneLengths.kb

for(i in colnames(rpkm)){
    message("Normalising for ", i)
    rpkm[,i] = rpkm[,i]/(dataNorm[i]*10^6)
}
rpkm_annot = mergeWithAnnotation(rpkm, annot_data) # merge rpkm=rownames and annot_data = gene_id
info(" Write RPPM to csv")
write.table(rpkm, file=paste(title,".rpkm.tsv",sep=""), sep="\t", quote=F)
write.csv(rpkm, file=paste(title,".rpkm.csv",sep=""))
write.csv(rpkm_annot, file=paste(title,"annot.rpkm.csv",sep=""), row.names=F)

save.image(paste0(title,"countTableToTPM.RData"))

#gene = "ENSG00000000938"
# > head(counttable)
# BCOR_shRNA1 BCOR_shRNA2 Control_shRNA1 Control_shRNA2 Control2_shRNA1 Control2_shRNA2 Ring1_RNF2_shRNA1
# ENSG00000223972.5           1           0              0              0               0               0                 0
# ENSG00000227232.5           1           0              0              0               1               0                 1
# ENSG00000278267.1           0           0              0              0               0               0                 0
# ENSG00000243485.5           0           0              0              0               2               0                 2
# ENSG00000284332.1           0           0              0              0               0               0                 0
# ENSG00000237613.2           0           0              0              0               0               0                 0
# Ring1_RNF2_shRNA2
# ENSG00000223972.5                 0
# ENSG00000227232.5                 0
# ENSG00000278267.1                 0
# ENSG00000243485.5                 1
# 
# > head(sample.matrix)
# condition    libtype
# BCOR_shRNA1     BCOR_shRNA paired-end
# BCOR_shRNA2     BCOR_shRNA paired-end
# Control_shRNA1     Control paired-end
# Control_shRNA2     Control paired-end
# Control2_shRNA1   Control2 paired-end
# Control2_shRNA2   Control2 paired-end
# 
# > head(annot_data)
# gene_name           gene_id           coords strand  source
# 1     DDX11L1 ENSG00000223972.5 chr1:11869-14409      +  HAVANA
# 2      WASH7P ENSG00000227232.5 chr1:14404-29570      -  HAVANA
# 3   MIR6859-1 ENSG00000278267.1 chr1:17369-17436      - ENSEMBL
# 4 MIR1302-2HG ENSG00000243485.5 chr1:29554-31109      +  HAVANA
# 5   MIR1302-2 ENSG00000284332.1 chr1:30366-30503      + ENSEMBL
# 6     FAM138A ENSG00000237613.2 chr1:34554-36081      -  HAVANA
# 
#res.results = getDESeq2(counttable, sample.matrix, annot_data, "Control", "BCOR_shRNA")


