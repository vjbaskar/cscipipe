# Desktop options
# title = "temp"
# sampleFile = "/Volumes/scratch119/ETIENNE/2018.03.Notch.Etienne.RNASeq/DATA/sampleTable.prj.tsv"
# cols = comments =  "SAMPLEID,SAMPLENAME"
# inputType="tsv"
# outputType="tsv"
# outpt="temp.tsv"
# source("../src/basic_functions.R")
#


#Given an tsv/csv file get the required columns

suppressPackageStartupMessages({
    library(dplyr)
})

mconf_installdir=Sys.getenv("mconf_installdir")
basic_functions_file = paste0(mconf_installdir,"/src/basic_functions.R")
source(basic_functions_file)


args_return=Sys.getenv("args_return")

helpVals = matrix(c(
	"t", "title", 1, "Title",
	"s","sampleFile", 1, "file to be filtered",
	"c","comments",1, "colnames to be filtered. eg [SAMPLEID,SAMPLENAME]",
	"0","cmd0",1, "Input file type [tsv,csv]",
	"1","cmd1", 1, "Output file type [tsv,csv]",
	"O","outpt",1, "Output file"
	),
	ncol = 4, byrow = T
)


if(args_return > 0){
    Help(helpVals)
    quit(status = 2)
}
Help(helpVals)

title = Sys.getenv("title")
sampleFile = Sys.getenv("sampleFile")
cols = Sys.getenv("comments")
inputType=Sys.getenv("cmd0")
outputType=Sys.getenv("cmd1")
outpt=Sys.getenv("outpt")

getPrjName <- function(projectCSV){
    temp = strsplit(projectCSV, "[.]")[[1]]
    temp = temp[-length(temp)]
    temp
    pname = paste0(temp, sep="", collapse = ".")
    return(pname)
}

message("INPUT PARAMS")
message("------------")
message("title = ", title)
message("samplefile = ", sampleFile)
message("cols = ", cols)
message("input type = ", inputType)
message("output type = ", outputType)
message("Output = ", outpt)

if(inputType=="tsv"){
	info("Reading file as tsv:")
	x = read.table(sampleFile, sep="\t", header=T, comment.char='#')
}
if(inputType=="csv"){
	info("Reading file as csv:")
	x = read.csv(sampleFile, header=T, comment.char='#')
}

x <- cleanSheet(x)
reqCols = strsplit(cols, ",")[[1]]
x = x[,reqCols]

if(outputType=="tsv"){
    info("Write out as tsv: ")
	x = write.table(x, file=outpt, sep="\t", quote=F ,row.names=F)
}
if(outputType=="csv"){
    info("Write out as csv: ")
	x = write.csv(x, file=outpt, row.names=F)
}