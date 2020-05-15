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
	"s","sampleFile", 1, "bam file list"
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


command=paste0(mconf_installdir,"/rmds/bamstats.Rmd")
args = commandArgs(trailingOnly=T)
p=getwd()
message(command)
params = list(fileName = sampleFile, wd=p)
params
rmarkdown::render(command, params = params, envir = new.env(parent = globalenv()), output_dir = paste0(p,"/bamstats/"), output_format = "all")