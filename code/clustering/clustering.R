#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]
output_file <- args[2]

set.seed(42)  

library(Biostrings)  
library(DECIPHER)    

dna <- readDNAStringSet(file_path)
dna

clusters <- Clusterize(myXStringSet = dna, cutoff = 0.1, processors = 8)

# Create data frame
clusters_df <- data.frame(
    Seq_ID = names(dna),
    Cluster_ID = clusters
)

# Rename the column
names(clusters_df)[2] <- "Cluster_ID"

write.csv(clusters_df, file = output_file, row.names = FALSE)
