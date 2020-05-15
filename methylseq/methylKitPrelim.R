library(methylKit)
f = "/lustre/scratch119/casm/team163gv/vm11/GONIA/2018.11.Gonia.DPMX_RRBS.MethylSeq/BISMARK/EXTRACTOR/DP_G7_OX.sormadup.bam"
args = commandArgs(trailingOnly = T)
organism = as.character(args[1])
fileName  = as.character(args[2])

f =  system(paste0("readlink -f ", fileName), intern = T)
f

getName <- function(f){
	temp = strsplit(f, "/")[[1]]
	gsub(".bam|.sam|.sormadup", "", temp[length(temp)], perl = T )
}

sample_id = getName(f)
sample_id

obj=processBismarkAln(f,  sample_id ,assembly=organism,save.folder=paste0(sample_id, "_pba"), save.context=c( "CpG","CHG","CHH"),read.context= c("CpG"), save.db = TRUE)

save.image(paste0(sample_id, ".RData"))

