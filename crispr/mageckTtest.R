# basic commands for an R program in cscipipe

library(getopt)
mconf_installdir=Sys.getenv("mconf_installdir")
source(paste0(mconf_installdir,"/src/basic_functions.R"))

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
            * Required Args:
            
            -t|--title
            -0|--cmd0: file for condn1 [eg. mageck_Plasmid_RapidRun..vs..NoDrug.sgrna_summary.txt]
            -1|--cmd1: file for condn2 [eg. mageck_Plasmid_RapidRun..vs..Drug.sgrna_summary.txt]
            -2|--cmd2: Pass cutoff for gRNA counts in at least one condn (C) 
            -3|--cmd3: Pass cutoff for gRNA FDR in at least one condn (P)
            -4|--cmd4: # of gRNA to pass both C and P to be considered for further analysis. 
                        If you want to run for only your genes of interest enter value -1.
            -f|--inputFile: list of genes to include
            ")
    quit(status = 2)
}

title = Sys.getenv("title")
c1_table = Sys.getenv("cmd0")
c2_table = Sys.getenv("cmd1")
countCutOff = Sys.getenv("cmd2")
fdrCutOff = Sys.getenv("cmd3")
sgRNACutoff = Sys.getenv("cmd4")
geneFile = Sys.getenv("inputFile")

# setwd("/Volumes/scratch117/ETIENNE/2018.05.PPM1D.EtienneJoanne.CRISPR/MAGECK")
# c1_table = "mageck_Plasmid_RapidRun..vs..C_2102NT.sgrna_summary.txt"
# c2_table = "mageck_Plasmid_RapidRun..vs..C_2102ARAMUT.sgrna_summary.txt"
# countCutOff = 100
# fdrCutOff = 0.05
# sgRNACutoff = 2
# geneFile = "temp.txt"

condn1_sgRNATable = read.table(c1_table, head=T)
condn2_sgRNATable = read.table(c2_table, head=T)


message("Condition 1: ", c1_table)
message("Condition 2: ", c2_table)
message("sgRNA count cutoff: ", countCutOff)
message("sgRNA FDR cutoff: ", fdrCutOff)
message("# sgRNA cutoff: ", sgRNACutoff)

condn1_sgRNATable$countMean <- paste0(round(condn1_sgRNATable$control_mean), "_", round(condn1_sgRNATable$treat_mean))
condn2_sgRNATable$countMean <- paste0(round(condn2_sgRNATable$control_mean), "_", round(condn2_sgRNATable$treat_mean))


clearFDR <- function(matx, fdrCutOff) {
    matx$goodFDR <- 0
    matx [ matx$FDR <= fdrCutOff, "goodFDR" ] <- 1
    return(matx)
}

clearCounts <- function(matx, countCutOff ) {
    matx$goodCount <- 0
    meancount = matx [ , c("control_mean", "treat_mean")]
    m <- apply(meancount, 1, max)
    matx [ m > countCutOff, "goodCount"] <- 1
    return(matx)
}


condn1_sgRNATable_sig = clearFDR(condn1_sgRNATable, fdrCutOff)
condn2_sgRNATable_sig = clearFDR(condn2_sgRNATable, fdrCutOff)

condn1_sgRNATable_sig = clearCounts(condn1_sgRNATable_sig, countCutOff)
condn2_sgRNATable_sig = clearCounts(condn2_sgRNATable_sig, countCutOff)

temp1 = unique(as.character(condn1_sgRNATable[ condn1_sgRNATable$FDR < fdrCutOff, "Gene"]))
temp2 = unique(as.character(condn2_sgRNATable[ condn2_sgRNATable$FDR < fdrCutOff, "Gene"]))

if(! geneFile == "" ){
    message("Given a gene file")
    userGenes = as.character(read.table(geneFile)$V1)
    genes = userGenes
    if(sgRNACutoff > -1){
        genes = unique(c(genes, temp1, temp2))    
    }
} else {
    message("No gene file")
    genes = unique(c(temp1, temp2))
}
message("Total genes to process = ", length(genes))

getGeneList <- function(condn1_table, condn2_table, genes){
    test_data = c()
    cat("Gene passed = ")
    allGenes = unique(condn1_table$Gene, condn2_table$Gene)
    genes = unique(allGenes [ allGenes %in% genes])
    for(geneName in genes){
        gene_condn1 <- condn1_table [ condn1_table$Gene == geneName, ]
        gene_condn2 <- condn2_table [ condn2_table$Gene == geneName, ]
        rownames(gene_condn2) <- gene_condn2$sgrna
        gene_condn2 <- gene_condn2 [ as.character(gene_condn1$sgrna), ]
        gene_goodFDR <- sum(gene_condn1$goodFDR + gene_condn2$goodFDR)
        gene_goodCount <- sum(gene_condn1$goodCount + gene_condn2$goodCount)
        if(gene_goodCount > sgRNACutoff & gene_goodFDR > sgRNACutoff ){
            cat(geneName," ")
            cond1 <- gene_condn1$LFC
            cond2 <- gene_condn2$LFC
            
            # Paired Student's t-test (Since the variance is not declared same it is Welch T-test). Two-sided
            gene_test = t.test(cond1, cond2, paired=T)
            temp = data.frame(gene = geneName, pval = gene_test$p.value, 
                              cond1_meanLFC = mean(round(cond1,2)),
                              cond2_meanLFC = mean(round(cond2,2)),
                              cond1_lfc = paste0(round(cond1,2),collapse="|"), 
                              cond2_lfc = paste0(round(cond2,2),collapse="|"), 
                              cond1_fdr = paste0(round(gene_condn1$FDR,3),collapse="|"), 
                              cond2_fdr = paste0(round(gene_condn2$FDR,3),collapse="|"), 
                              cond1_counts =  paste0(gene_condn1$countMean,collapse="|"), 
                              cond2_counts =  paste0(gene_condn2$countMean,collapse="|")
                              
            )
            test_data = rbind(test_data, temp)
        } 
    }
    return(test_data)
}
test_data = c()
if(sgRNACutoff > -1){ 
    test_data <- getGeneList(condn1_sgRNATable_sig, condn2_sgRNATable_sig, genes)
}
if(! geneFile == ""){
    message("Running the given gene list separately")
    condn1_sgRNATable$goodFDR <- 1 ; condn1_sgRNATable$goodCount <- 1
    condn2_sgRNATable$goodFDR <- 1 ; condn2_sgRNATable$goodCount <- 1
    temp = getGeneList(condn1_sgRNATable, condn2_sgRNATable, userGenes)
}
test_data = rbind(test_data, temp)

test_data$FDR <- p.adjust(test_data$pval)
#
write.table(test_data, paste0(title, ".mageckTtest.tsv"), sep="\t", quote=F, row.names=F)
save.image(paste0(title, ".mageckTtest.RData"))
cat("\n")
# # Testing
# for(i in  seq(-5,5)){
#     cond1 = rnorm(5, mean = 2)
#     cond2 = rnorm(5, mean = i)
#     message(i)
#     t <- t.test(cond1, cond2, paired=T)
#     print(t$p.value)
#    w <-  wilcox.test(cond1, cond2, paired=T)
#    print(w$p.value)
#     
# }

