
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
	"s", "sampleFile", 1, "Input fastq file",
	"0","cmd0",1, "Library file. Look for it in ..../SHARED/reference/crispr/libs/"
	),
	ncol = 4, byrow = T
)

if(args_return > 0){
    Help(helpVals)
    quit(status = 2)
}
Help(helpVals)


title = Sys.getenv("title")
fqFile = Sys.getenv("sampleFile")
lib=Sys.getenv("cmd0")



suppressPackageStartupMessages({
    library(ShortRead)
    library(pryr)
    library(ggplot2)
    library(ggthemes)
    library(ggpubr)
    library(gridExtra)
})

count_na = 25
f=gsub(".fq+.gz","", fqFile)
ofile <- paste0(f, ".counts.tsv")

cat("Reading library ... ", lib, "\n")
libData = read.table(lib, head=T)
libData = dplyr::as_tibble(libData)


cat("Streaming in fq file ", fqFile)
countsData <- c()
fqStr <- ShortRead::FastqStreamer(fqFile)
i=0
repeat {
    fq <- ShortRead::yield(fqStr)
    if (length(fq) == 0)
        break
    i = i+1
    message("Read ", i, 'x', length(fq), " reads", ". Memory used: ", round(mem_used()/10^6), " MB")
    temp = dplyr::as_tibble(as.data.frame(table(sread(fq))))
    countsData = bind_rows(countsData, temp)
}


"Counting gRNA ..."
colnames(countsData) <- c("gRNA_ID", "n")
mem_change(count_data <- group_by(countsData, gRNA_ID) %>% mutate(count = sum(n))) 
count_data <- select(count_data, -n)
count_data <- unique(count_data)
mem_used()

# " gRNA counting ..."
# creads = do.call("c", unlist(creads))
# system.time(count_data <-  table(creads))
# count_data = as_tibble(count_data)
# colnames(count_data) <- c("gRNA_ID","count")

"Merging gRNA counts with library"
merged_counts <- full_join(libData, count_data, by = c("Guide_sequence" = "gRNA_ID"))
merged_counts <- filter(merged_counts, (is.na(gRNA_ID) & count > count_na) | (! is.na(gRNA_ID)))
" -- if gRNA count is NA then substitute with 0"
merged_counts <- mutate(merged_counts, count = ifelse(is.na(count), 0, count))
mem_used()

"Clean up merged data"
merged_counts = merged_counts %>% mutate(Gene = case_when( is.na(Gene) ~ "UN", TRUE ~ as.character(Gene) ), gRNA_ID = case_when( is.na(gRNA_ID) ~ "UN", TRUE ~ as.character(gRNA_ID) ))



"Plot some preliminary stats"
merged_counts_stats = merged_counts %>% mutate(Annotation = case_when(Gene == "UN" ~ "unknown", Gene != "UN" ~ "known"))
count_stats <- merged_counts_stats %>% group_by(Annotation) %>% summarise(counts = sum(count))
complexity  <- data.frame(rank = 1:nrow(merged_counts), log2counts = sort(log2(merged_counts$count+1)))
top_unk <- head(merged_counts_stats %>% filter(Gene == "UN") %>% arrange(desc(count) ), 10)

pdf(paste0(f, ".crispr.pdf"))
p_complexity <- ggplot(data = complexity ) + aes(x = rank, y = log2counts) + geom_line(colour = "red", size=2) + xlab("gRNA Rank") + ylab("log2 (# reads + 1)") + theme_igray()+ ggtitle(f)
p_counts <- ggplot(data = count_stats) + aes(x = Annotation, y = counts, fill = Annotation) + geom_bar(stat="identity") + scale_fill_gdocs() + theme_excel_new() + xlab("gRNA") + ylab("# Reads") 
p_table <- ggtexttable(top_unk[,c("Guide_sequence", "count")],theme = ttheme("mOrangeWhite"), rows = NULL, cols = c("gRNA unk", "n"))
ggarrange(p_complexity, ggarrange(p_counts,p_table , ncol = 2), nrow=2)
dev.off()


"Write out data"
colnames(merged_counts)[ncol(merged_counts)] <- f
write.table(merged_counts, sep="\t", row.names=F, quote	=F, file = ofile)
save.image( file = paste0(f,".counts.RData"))


