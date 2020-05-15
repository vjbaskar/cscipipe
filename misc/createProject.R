#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library( getopt )
})


spec = matrix(c(
    "sampleSheet" ,  "s" , 1, "character",
    "cols", "c", 1, "character",
    "tsv", "t" , 1, "character",
    "help",  "h" , 0, "logical"
    ), byrow=TRUE, ncol=4);

opt = getopt(spec)
if ( !is.null(opt$help)  ) {
    message("\n")
    message(getopt(spec, usage=TRUE),"\n");
    q(status=1);
}

if ( length(opt) == 1  ) {
    message("\n")
    message(getopt(spec, usage=TRUE),"\n");
    q(status=1);
}

projectCSV = opt$sampleSheet
message("Prj table = ", projectCSV)
cnFields = opt$cols
message("Cols requested  = ", cnFields)
cnFields = strsplit(cnFields, ",")[[1]]
tsv = opt$tsv
message("File is tsv  = ", tsv)


#projectCSV = "2017.10.K562_Sf3b1Mutant_crispr.csv"
pwd = getwd()



getPrjName <- function(projectCSV){
    temp = strsplit(projectCSV, "[.]")[[1]]
    temp = temp[-length(temp)]
    temp
    pname = paste0(temp, sep="", collapse = ".")
    return(pname)
}

info <- function(text){
    text =  sapply(text, function(x) paste0(x, sep=" ", collapse = " "))
    d = date()
    cat(" [ Info: ", d, "]\n")
    cat(text,"\n")
}

errorMessage <- function(text){
    text =  sapply(text, function(x) paste0(x, sep=" ", collapse = " "))
    d = date()
    cat(" [ Error: ", d, "]\n")
    cat(text,"\n")
    cat("Quitting")
    message(" [ Error: ", d, "]\n")
    message(text,"\n")
    cat("Quitting")
    quit(save="yes", status=1)
}

createFileLink <- function(fileName){
    return(paste0(projectDir,"/", fileName))
}

stripNALines <- function(ss){
    # row wise. If Run
    temp = apply(ss, 2, unique)
    ss = ss[, - which(is.na(temp))]
    ss = ss [ !is.na(ss[,2]), ]
    ss = ss [ !is.na(ss[,1]), ]
    
    return(ss)
}

cleanSheet <- function(ss){
    ss = stripNALines(ss)
    for(j in 1:ncol(ss)){
        ss[,j] <- as.character(ss[,j])
    }
    
    for(i in 1:nrow(ss)){
        for(j in 1:ncol(ss)){
            t <- as.character(ss[i,j])
            keepOnlyAlphaNumeric(t)
            
            ss[i,j] <- keepOnlyAlphaNumeric(t)
        }
    }
    cn <- colnames(ss)
    cn <- sapply(cn, keepOnlyAlphaNumeric)
    cn <- toupper(cn)
    colnames(ss) <- cn
    return(ss)
}



keepOnlyAlphaNumeric <- function(text){
    specialChars <- "a1~!@#$%^&*(){}_+:\"<>?./;'[]-= "
    # text = gsub("[[ ].$]", "",text)
    text = gsub("[ ]", "_",text)
    text = gsub("[^[:alnum:]^,^_^-]", "",text)
    text = gsub("_*$", "",text)
    return(text)
}

projectName = getPrjName(projectCSV)

#info(c("Project Name=", projectName))

#projectDir = paste0(pwd, "/", projectName)
#info(c("Project Directory is = ", projectDir))

# if(dir.exists(projectDir)) { 
#     errorMessage(c("Project Directory exists")) 
# } else {
#     dir.create(projectDir)
# }

if(tsv==1) { 
	sampleSheet = read.table(projectCSV, blank.lines.skip = TRUE, comment.char = "#", fill = NA)
	} else { 
	sampleSheet = read.csv(projectCSV, blank.lines.skip = TRUE, comment.char = "#", fill = NA)
}

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




