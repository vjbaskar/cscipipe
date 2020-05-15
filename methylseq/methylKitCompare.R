suppressPackageStartupMessages({
    library(getopt)
})
outpt_prefix = "methylkit"
#outpt_rdata = paste0(outpt_prefix, ".RData")
#outpt_tiles = paste0(outpt_prefix, ".tiles.tsv")
#outpt_res = paste0(outpt_prefix, ".res.tsv")
sampleFile = "samples.test.ox.txt"
qcut = 0.05
percCutOff = 10
organism = "mm10"
ncores = 1
lowcountCutoff = 10
tilingWindowSize = 200
tilingStepSize = 50
#cpgBedFile ="../ANNOT/mm10.cpg.bed"

spec = matrix(c( 
    'outpt_prefix', 'O', 1, "character", 
    'sampleFile' , 'f', 1, "character", 
    'qcut' , 'v', 1, "double", 
    'percCutOff' , 'c', 1, "double", 
    'lowcountCutoff','l',1, "integer",
    'tilingWindowSize','w',1,"integer",
    'tilingStepSize','s',1, "integer",
    'organism' , 'o', 1, "character",
    'cpgBedFile','b',1,"character",
    'ncores' , 'p', 1, "integer",
    'help', 'h', 0, "integer"
), byrow=TRUE, ncol=4); 
opt = getopt(spec)

if ( !is.null(opt$help) ) { 
    cat(getopt(spec, usage=TRUE)); 
    q(status=1); 
}

suppressPackageStartupMessages({
    library(methylKit)
    library(genomation)
    library(dplyr)
    library(annotatr)
})

outpt_prefix = opt$outpt_prefix
outpt_rdata = paste0(outpt_prefix, ".RData")
outpt_tiles = paste0(outpt_prefix, ".tiles.tsv")
outpt_res = paste0(outpt_prefix, ".res.tsv")
outpt_cpg = paste0(outpt_prefix, ".cpg.tsv")
outpt_pdf = paste0(outpt_prefix, ".mkit.pdf")

sampleFile = opt$sampleFile
qcut = opt$qcut
percCutOff = opt$percCutOff
organism = opt$organism
ncores = opt$ncores
lowcountCutoff = opt$lowcountCutoff
tilingWindowSize = opt$tilingWindowSize
tilingStepSize = opt$tilingStepSize
cpgBedFile = opt$cpgBedFile

# Read sample file
message("Reading sample file = ", sampleFile)
x = read.table(sampleFile)
colnames(x) <- c("files", "names", "treatment")
file.list=as.list(as.character(x$files))
sampleids = as.list(as.character(x$names))
t = as.list(x$treatment)

# Create methylkit objects
message("Creating methyl raw list objects")
methRawList_data=methRead(file.list,
                          sample.id=sampleids,
                          assembly=organism,
                          treatment=as.numeric(x$treatment),
                          context="CpG"
)


# methTabixList_data=methRead(file.list,
#               sample.id=sampleids,
#               assembly=organism,
#               treatment=t,
#               context="CpG",
# 			  dbtype = "tabix",
#               dbdir = "methylDB"
#               )

save.image(outpt_rdata)
##### Common Functions

od="MN"
stest="F"
padj = "BH"
mc.cores=ncores
calcDiff <- function(x){
    temp = calculateDiffMeth(x, overdispersion=od, test=stest, adjust = padj, mc.cores=ncores)
    return(temp)
}

lowcountCutoff = lowcountCutoff
lo.perc=NULL
hi.count=NULL
hi.perc=99.9
filtByCov <- function(mRawList){
    filterByCoverage(methRawList_data,lo.count=lowcountCutoff,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
    
}


getNormCounts <- function(mRawList, diffdata){ # mRawList = cpg_counts # diffdata = diffmeth_cpg
    temp1 <- as(diffdata, "data.frame")
    counts_norm = methylKit::normalizeCoverage(mRawList)
    unite_cpg_norm = methylKit::unite(counts_norm)
    temp = as(unite_cpg_norm, "data.frame")
    sample_ids = unite_cpg@sample.ids 
    cn <- lapply(sample_ids, function(x) paste0(x, c("_cov","_Cs","_Ts")))
    cn <- unlist(cn)
    colnames(temp)
    cn <- c(c("chr", "start","end","strand"), cn)
    colnames(temp) <- cn
    temp <-left_join(temp1, temp)
    return(temp)
}


message("#### Running per base level")
# Filter by sample cov
message("Filter coverage by low counts = ", lowcountCutoff)
meth_filt=filtByCov(methRawList_data)
# Unite data
message("Unite methykit data")
meth_unite=unite(meth_filt, destrand=FALSE)
# Unite by groups
message("Unite methylkit data by conditions")
meth_condn = unite(meth_filt, min.per.group=1L)

message("#### Running per tile level")
message("Tiling methyl counts")
tiles=tileMethylCounts(meth_filt,win.size=tilingWindowSize,step.size=tilingStepSize)
#head(tiles[[1]],3)
message("Filtering tiles")
#tiles_filt = filtByCov(tiles)
tiles_filt = tiles
tiles_unite = unite(tiles_filt, destrand = FALSE)

message("Saving data ...")
save.image(outpt_rdata)

#message("#### Running CpG islands")

#message("Reading CpG file ...")
#cpg.obj = genomation::readFeatureFlank(cpgBedFile,feature.flank.name=c("CpGi","shores"))

#message("Filtering and tiling CpG islands ...")
#meth_filt_cpg = methylKit::selectByOverlap(meth_filt, cpg.obj$CpGi)
#tiles_cpg = tileMethylCounts(meth_filt_cpg,win.size=tilingWindowSize,step.size=tilingStepSize)

#message("Getting region counts per CpG ...")
#cpg_counts = regionCounts(meth_filt, cpg.obj$CpGi)
#message("Uniting these ...")
#unite_cpg = unite(cpg_counts)




# diff meth
message("Calc diff meth ...")
message("Calc diff meth in bp level")
diffmeth = calcDiff(meth_unite)
message("Calc diff meth in tiles level")
diffmeth_tiles = calcDiff(tiles_unite)
#message("Calc diff meth at CpG level")
#diffmeth_cpg = calcDiff(unite_cpg)

# Get diffmeth data

message("Get diff meth data")
diffmeth_res = getMethylDiff(diffmeth,difference=percCutOff,qvalue=qcut) # percCutOff= 0-100 ; qcut = 0 - 1; type = hyper|hypo|all 
diffmeth_tiles = getMethylDiff(diffmeth_tiles,difference=percCutOff,qvalue=qcut)
#diffmeth_cpg = getMethylDiff(diffmeth_cpg,difference=percCutOff,qvalue=qcut)




message("Saving data ...")
save.image(outpt_rdata)


message("Creating comprehensive output")

#annotations = build_annotations(genome = organism, annotations = paste0(organism,"_basicgenes"))
#annotations = build_annotations(genome = organism, annotations = paste0(organism,"_cpgs"))

getAnnotationData <- function(countsdata, diffmethdata){ # countsdata = cpg_counts ; # dimethdata = diffmeth_cpg
    temp <- getNormCounts(countsdata, diffmethdata)
    x = as(temp, "GRanges")
    mapping = annotate_regions(x, annotations)
    df_mapping = data.frame(mapping)
    return(df_mapping)
}




message("Write output")
write.table(diffmeth_res, file = outpt_res, sep="\t", quote=F, row.names = F)

temp <- getAnnotationData(tiles_filt, diffmeth_tiles)
write.table(temp, file = outpt_tiles, sep="\t", quote=F, row.names = F)
#pdf(outpt_pdf,7,8)
#par(mar = c(5,5,5,5))
#pie(table(temp$annot.type), main = "Sliding Window Mode")
#dev.off()
#temp = getAnnotationData(cpg_counts, diffmeth_cpg)
#write.table(diffmeth_cpg, file = outpt_cpg, sep="\t", quote=F, row.names = F)

# Correlation, plotting and pca
#message("PCA")
#pdf("corr.pdf")
#getCorrelation(meth_unite, plot=T)
#clusterSamples(meth_unite, dist="correlation", method="ward", plot=TRUE)
#dev.off()

#hc = clusterSamples(meth_unite, dist="correlation", method="ward", plot=FALSE)

#PCASamples(meth_unite, screeplot=TRUE)

#PCASamples(meth_unite)

