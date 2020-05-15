suppressPackageStartupMessages({
    library(getopt)
})


spec = matrix(c( 
    'rdatafile' , 'f', 1, "character", 
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

# Read data file
load(rdatafile)

