# basic commands for an R program in cscipipe
#### In DiffBind the logFC = A/B not B/A

suppressPackageStartupMessages({
    library(getopt)
    library(dplyr)
    library("DiffBind")
})


title = Sys.getenv("title")
filename = Sys.getenv("inputFile")
genome = Sys.getenv("organism")
g1 = Sys.getenv("cond1")
g2 = Sys.getenv("cond2")
ofold=Sys.getenv("outpt")
ofold=paste0(ofold,"/",g1,"_vs_",g2)

# setwd("~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_ChIPSeq/DIFFBIND/H3K4ME3")
# filename = "samples.diffbind.txt"
# genome = "mm10"
# g1 = "gRNA_empty"
# g2 = "gRNA_Setd1b"


message("Loading organism data ")
#### Organism data
if(genome == "hg38"){
    message("Loading packages for hg38 ....")
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    annoDb="org.Hs.eg.db"
    
}

if(genome == "hg19"){
    message("Loading packages for hg19 ....")
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    annoDb="org.Hs.eg.db"
    
}

if(genome == "mm10"){
	message("Loading packages for hg19 ....")
	library(org.Mm.eg.db)
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
	annoDb="org.Mm.eg.db"

}

message("==== Conditions ====")
g1
g2

samples = read.table(filename, head=T, sep="\t")
message("==== Sample table ====")

samples = samples %>% filter(Condition %in% c(g1,g2))
samples

#ofold = gsub(".txt", "outpt", filename)
system(paste0("mkdir ", ofold), intern = F)

# Filter
fdrcutoff = 0.10
fc_up = 0.5
fc_down = -0.5




# Read peaksets
message("Reading peaks ...")
peakSets = dba(sampleSheet = samples)
peakSets

# Affinity scores!
message("Getting counts")
#counts = dba.count(peakSets, bRemoveDuplicates = T, bParallel = 2)
counts = dba.count(peakSets, bRemoveDuplicates = T)
counts

# Set contrasts
counts.bkup = counts
#counts <- dba.contrast(counts, group1 = counts$masks$dmso, group2 = counts$masks$ibet, name1 = "dmso", name2 = "ibet")
counts <- dba.contrast(counts, group1 = counts$masks[[g1]], group2 = counts$masks[[g2]], name1 = g1, name2 = g2)

# Diff binding testing

diffbind_stats = dba.analyze(counts, method=DBA_DESEQ2)


diffbind_stats
diffbind_db = dba.report(diffbind_stats, th=1, bCounts=TRUE, bCalled=TRUE, file="diffbind")


# Plots 
pdf(paste0(ofold,"/diffbind.stats.pdf"))
plot(diffbind_stats, main = "Normalised count correlation")
plot(peakSets, main = "peak correlation")
dev.off()
# FDR cutoff
filt_data = diffbind_db [ diffbind_db$FDR <= fdrcutoff ]

# Down 
filt_data_down = filt_data [ filt_data$Fold <= fc_down ]
length(filt_data_down)

# Up 
filt_data_up = filt_data [ filt_data$Fold >= fc_up ]
length(filt_data_up)


# Map peaks to genes

library(ChIPseeker)

mapPeaks <- function(x, txdb = txdb, annoDb = annoDb){
    peakAnno <- annotatePeak(x, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb = annoDb , overlap = "all")
    cn <- colnames(mcols(x))
    gr = as.GRanges(peakAnno)
    cn <- c(cn, c("SYMBOL","distanceToTSS","annotation", "GENENAME"))
    gr = gr[, cn]
    #names(diffbind_db_annot) <- paste0("diffPeaks_",1:length(diffbind_db_annot))
    return(list=c(peakAnno,gr))
}



map.filt_data = mapPeaks(filt_data, txdb, annoDb)
map.filt_data_down = mapPeaks(filt_data_down, txdb, annoDb)
map.filt_data_up = mapPeaks(filt_data_up, txdb, annoDb)
map.all =  mapPeaks(diffbind_db, txdb, annoDb)


plotDistn <- function(x, pdffile){
    library(ggplot2)
    library(ggthemes)
    pdf(pdffile,9, 2.7)
    p <- plotAnnoBar(x)
    p <- p + scale_fill_calc() + theme(plot.background = element_rect(fill = 'white'), panel.background  = element_rect(fill = 'white'))
    print(p)
    p <- plotDistToTSS(x, title="Distribution of peaks relative to TSS", cols = rainbow(9))
    p <- p + scale_fill_calc() + theme_dark() + theme(plot.background = element_rect(fill = 'white'), panel.background  = element_rect(fill = 'white'))
    print(p)
    dev.off()
}


plotDistn(map.filt_data_down[[1]], paste0(ofold,"/downreg.pdf"))
plotDistn(map.filt_data_up[[1]], paste0(ofold,"/upreg.pdf"))
plotDistn(map.filt_data[[1]], paste0(ofold,"/both.pdf"))

logfile=paste0(ofold,"/log.txt")
cat("FDR cutoff = ", fdrcutoff, "\n", file = logfile, append=F)
cat("FC cutoff = ", abs(fc_down), "\n", file = logfile, append=T)
cat("Total regions = ", length(map.filt_data[[2]]), "\n", file = logfile, append=T)
cat("Total up regions = ", length(map.filt_data_up[[2]]), "\n", file = logfile, append=T)
cat("Total down regions = ", length(map.filt_data_down[[2]]), "\n", file = logfile, append=T)
cat("diffReg table = ", filename, "\n", file = logfile, append=T)
write.table(samples, file = logfile, append=T)


# Final data
finalData = as.data.frame(map.all[[1]])
finalData$status <- "nochange"
finalData[ finalData$FDR <= fdrcutoff & finalData$Fold <= fc_down, "status"] <- "down"
finalData[ finalData$FDR <= fdrcutoff & finalData$Fold >= fc_up, "status"] <- "up"
write.csv(finalData, file=paste0(ofold,"/diffbind.csv"), row.names=F)

# Write bed
finalData$name = paste0("diffPeak.", finalData$SYMBOL, ":", finalData$distanceToTSS, "..FC-FDR:", finalData$Fold, "-",finalData$FDR )
finalData$log10Q = -log10(finalData$FDR)

writeBED <- function(x, fname){
     reqCols = c("seqnames", "start", "end", "name", "log10Q","strand")
    x = x[,reqCols]
    write.table(x, file=paste0(ofold,"/",fname), sep="\t", quote=F, row.names=F, col.names=F)
}

writeBED(finalData[finalData$status == "down", ], "downreg.bed")
writeBED(finalData[finalData$status == "up", ], "upreg.bed")
writeBED(finalData[finalData$status == "nochange", ], "nochange.bed")




save.image(paste0(ofold,"/diffbind.RData"))
