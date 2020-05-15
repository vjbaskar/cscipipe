# basic commands for an R program in cscipipe

library(getopt)
mconf_installdir=Sys.getenv("mconf_installdir")
source(paste0(mconf_installdir,"/src/basic_functions.R"))

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
* Required Args:

-t|--title
-s|--sampleFile: master sample file csv
-c|--comments: colnames eg. SAMPLEID,SAMPLENAME

		")
    quit(status = 2)
}

title = Sys.getenv("title")
sampleFile = Sys.getenv("sampleFile")
cols = Sys.getenv("comments")



projectCSV = sampleFile
message("Prj table = ", projectCSV)
cnFields = cols
message("Cols requested  = ", cnFields)
cnFields = strsplit(cnFields, ",")[[1]]


#projectCSV = "2017.10.K562_Sf3b1Mutant_crispr.csv"
pwd = getwd()

projectName = getPrjName(projectCSV )

#info(c("Project Name=", projectName))

#projectDir = paste0(pwd, "/", projectName)
#info(c("Project Directory is = ", projectDir))

# if(dir.exists(projectDir)) { 
#     errorMessage(c("Project Directory exists")) 
# } else {
#     dir.create(projectDir)
# }

sampleSheet = read.csv(projectCSV, blank.lines.skip = TRUE, comment.char = "#", fill = NA)
sampleSheet <- cleanSheet(sampleSheet)
write.table(sampleSheet, file=paste0(projectName, ".prj.tsv"), sep="\t", quote=F, row.names=F)

#cnFields = c("RunID", "SampleID")
cnFields = toupper(cnFields)

checkFields <- function(ss, cn) {
    ss.cn <- colnames(ss)
    temp <- cn [ ! cn %in% ss.cn ]
    if(length(temp) > 0 ) { 
        message("[ Error ] Cannot find columns you have requested ....", temp)
        message("Aborting. Error val = 1")
        quit(status = 1)
        }
}
checkFields(sampleSheet, cnFields)
samples.extract = as.data.frame(sampleSheet[,cnFields])
print(samples.extract, row.names=F)
