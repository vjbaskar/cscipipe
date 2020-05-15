# basic commands for an R program in cscipipe

suppressPackageStartupMessages({
    library(getopt)
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