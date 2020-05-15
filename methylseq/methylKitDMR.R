#load("oxBS3.RData")
#load("BS2.RData")
#mc.cores = 4
#qcut = 0.05
#percCutOff = 10
#organism = "mm10"
#ncores = 4
#lowcountCutoff = 10
#tilingWindowSize = 200
#tilingStepSize = 50

mconf_installdir=Sys.getenv("mconf_installdir")
basic_functions_file = paste0(mconf_installdir,"/src/basic_functions.R")
source(basic_functions_file)

suppressPackageStartupMessages({
    library(getopt)
    
})

spec = matrix(c(
	'input_RData', 'i', 1, "character", 
	'bpWidth','w', 1, 'integer',
	'help', 'h', 0, "integer"
), byrow=TRUE, ncol=4); 
opt = getopt(spec)

if ( !is.null(opt$help) ) { 
    cat(getopt(spec, usage=TRUE)); 
    q(status=1); 
}
input_RData = opt$input_RData
bpWidth = opt$bpWidth
load(input_RData)

outpt_dmr = paste0(outpt_prefix, ".DMR.tsv")
outpt_dmr_bed = paste0(outpt_prefix, ".DMR.bed")
outpt_dmrAll_bed = paste0(outpt_prefix, ".DMR_all.bed")


info("Loading packages")

suppressPackageStartupMessages({
    library(methylKit)
    library(genomation)
    library(dplyr)
    library(annotatr)
    library(ChIPseeker)
    library(rtracklayer)
    library(ChIPseeker)
 
})

filtByCov <- function(mRawList, ...){
   temp =  filterByCoverage(methRawList_data,lo.count=lowcountCutoff,lo.perc=NULL, hi.count=NULL,hi.perc=100, ...)
    return(temp)
}


calcDiff <- function(x, ...){
    temp = calculateDiffMeth(x, overdispersion=od, test=stest, adjust = padj, mc.cores=ncores, ...)
    return(temp)
}


#####
#### Base level computation
# Filter count bases
info("Filter by coverage")
info(paste0("count cut off = ", lowcountCutoff))

methBP_filt = filtByCov(methRawList_data)
info("Raw count/coord = ")
dim(methRawList_data[[1]])
info("After filtering count/coord = ")
dim(methBP_filt[[1]])

# Unite the counts into a single obj

methBP_unite = unite(methBP_filt, destrand = F)
dim(methBP_unite)


# Convert to GRanges
info("Merging: bp level to region level")
methBP_gr = as(methBP_unite, "GRanges")
methBP_gr <- sortSeqlevels(methBP_gr)
methBP_gr <- sort(methBP_gr)

# If two good count bases are within 100bp of each other. Then merge them
temp = GenomicRanges::resize(methBP_gr, width = bpWidth, fix = "center")
methBP_gr_red = unique(reduce(temp, ignore.strand=T))
cat("Total DMRs shortlisted for counting = ", length(methBP_gr_red))

#### Region level computation
# counting in defined regions
info("Counting for DMRs")
methDMR = regionCounts(methBP_unite,methBP_gr_red)

info("Calc diff meth at DMR level using F test")
od="MN"
stest="F"
padj = "BH"
mc.cores=ncores
diffDMR_Fstat = calcDiff(methDMR)

info("Calc diff meth at DMR level using Chi square test (recommended)")
od="MN"
stest="Chisq"
padj = "BH"
mc.cores=ncores
diffDMR_Chisq = calcDiff(methDMR)
diffDMR_Chisq_gr = as(diffDMR_Chisq,"GRanges")
#write.table(as.data.frame(diffDMR_Chisq_gr), paste0(outpt_prefix, ".DMR.tsv"), row.names = F, quote=F, sep="\t", col.names = T)


#### Writing output
# write.table(as.data.frame(diffDMR_Chisq_gr), paste0(outpt_prefix, ".DMR.tsv"), row.names = F, quote=F, sep="\t", col.names = T)

# diffDMR_gr = as(diffDMR, "GRanges")
# names(diffDMR_gr) = paste0(paste0("dmr",sprintf("%07.0f",seq(length(diffDMR_gr)))))
# #diffDMR_gr$name = names(diffDMR_gr)
# write.table(as.data.frame(diffDMR_gr), outpt_dmr, row.names = F, quote=F, sep="\t", col.names = F)

info("Mapping to the genes")

if(organism == "mm10") { 
   	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
	annoDb = "org.Mm.eg.db"
}
if(organism == "hg38") { 
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	annoDb = "org.Hs.eg.db"
}

peakAnno <- annotatePeak(diffDMR_Chisq_gr, TxDb = txdb, annoDb = annoDb)
peakAnno_df = as.data.frame(peakAnno)
peakAnno_df %>% filter(SYMBOL == "Hoxa5")

write.table(peakAnno_df, outpt_dmr, row.names = F, sep="\t", quote = F)
save.image(paste0(outpt_prefix,".DMR.RData"))


peak_bed <- peakAnno_df %>% 
    mutate(chr = seqnames, start = start, end = end, name = paste0(outpt_prefix, seq(nrow(peakAnno_df))), score = -log10(qvalue), strand = ".") %>% 
    filter(qvalue <= qcut) %>% 
    dplyr::select(chr, start, end, name, score, strand)

writeLines(paste0("track type=bed name=",outpt_prefix," description='",outpt_prefix," DMR'"), con = outpt_dmr_bed)
write.table(peak_bed, sep="\t", quote=F, row.names = F, col.names = F, append = T, file = outpt_dmr_bed)

peak_bed <- peakAnno_df %>% 
    mutate(chr = seqnames, start = start, end = end, name = paste0(outpt_prefix, seq(nrow(peakAnno_df))), score = -log10(qvalue), strand = ".") %>% 
    dplyr::select(chr, start, end, name, score, strand)
writeLines(paste0("track type=bed name=all_",outpt_prefix," description='",outpt_prefix," DMR All with No qcut'"), con = outpt_dmrAll_bed)
write.table(peak_bed, sep="\t", quote=F, row.names = F, col.names = F, append = T, file = outpt_dmrAll_bed)


# #### Tiling data
# 
# meth_tile = tileMethylCounts(methRawList_data, win.size=tilingWindowSize,step.size=tilingStepSize,  mc.cores = mc.cores)
# meth_unite = unite(meth_tile, destrand = T)
# meth_unite_df = as(meth_unite,"data.frame") 
# meth_unite_gr = as(meth_unite,"GRanges") 
# 
# save.image()

# chr6:52,203,672-52,204,123
# chr6:52,203,557-52,204,187chr6:52,201,865-52,205,941

#hoxa5 = GRanges(seqnames = "chr6", ranges = IRanges(52201865, 52205941))
#temp = subsetByOverlaps(methBP_gr_red, hoxa5)
#temp

#temp = subsetByOverlaps(methBP_gr_red, hoxa5)
#temp


#x = as(methDMR, "GRanges")
#temp = subsetByOverlaps(x, hoxa5)
#temp

#temp = subsetByOverlaps(diffDMR1_gr, hoxa5)
#temp

