library( "DESeq2" )
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
first_arg <- args[1]  # This will access 'arg1'

PATH_COUNT="genecounts.csv"
PATH_META="metadata.csv"
PATH_OUTPUT_NORM="norm.csv"

setwd(first_arg)
countData <- read.csv(PATH_COUNT, header = TRUE, sep = ",")
head(countData)



metaData <- read.csv(PATH_META, header = TRUE, sep = ",")


dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy = TRUE)
dds <- DESeq(dds)


normalized_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(normalized_counts), 
          file=PATH_OUTPUT_NORM)