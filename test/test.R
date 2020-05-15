#!/usr/bin/env Rscript



library(getopt)

args_return=Sys.getenv("args_return")
if(args_return > 0){
	message("Help: ")
	quit(status = 2)
}

title = Sys.getenv("title")
message(title)



