set.seed(42)  

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")

library(Biostrings)  # For reading and handling DNA sequences
library(DECIPHER)    # For clustering sequences

# Specify the path to your FASTA file
file_path <- "genes.fasta"

# Read the DNA sequences from the file
dna <- readDNAStringSet(file_path)

# View the first few sequences to check if they were loaded correctly
dna

# Cluster the sequences
clusters <- Clusterize(myXStringSet = dna, cutoff = 0.1, processors = 8)

# Save the clusters to a CSV file
write.csv(clusters, file = "clusters_cutoff_0.1.csv")
