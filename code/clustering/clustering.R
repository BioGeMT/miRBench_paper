#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]
output_file <- args[2]
cutoff <- if (length(args) >= 3) as.numeric(args[3]) else 0.1
processors <- if (length(args) >= 4) as.integer(args[4]) else 8

if (is.na(cutoff)) {
    stop("Cutoff must be a numeric value.")
}

if (is.na(processors) || processors < 1) {
    stop("Processors must be a positive integer.")
}

set.seed(42)  

library(Biostrings)  
library(DECIPHER)    

dna <- readDNAStringSet(file_path)
dna

clusters <- Clusterize(myXStringSet = dna, cutoff = cutoff, processors = processors)

# Create data frame
clusters_df <- data.frame(
    Gene_ID = names(dna),
    Cluster_ID = clusters
)

write.csv(clusters_df, file = output_file, row.names = FALSE)
