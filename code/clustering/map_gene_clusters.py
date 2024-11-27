import argparse
import pandas as pd
from Bio import SeqIO

def detect_sequence_column(df, sequences):
    for col in df.columns:
        if any(seq in df[col].values for seq in sequences):
            return col
    raise ValueError("No column containing sequences found in gene file")

def load_clusters(cluster_csv):
    clusters_df = pd.read_csv(cluster_csv)
    
    if len(clusters_df.columns) == 2:
        if clusters_df.columns[0] == "" and clusters_df.columns[1] == "cluster":
            clusters_df.columns = ["GeneID", "ClusterID"]
        else:
            clusters_df = clusters_df[["Unnamed: 0", clusters_df.columns[-1]]]
            clusters_df.columns = ["GeneID", "ClusterID"]
    else:
        cluster_data = clusters_df[["Unnamed: 0", clusters_df.columns[-1]]]
        clusters_df = cluster_data
        clusters_df.columns = ["GeneID", "ClusterID"]
    
    return dict(zip(clusters_df["GeneID"], clusters_df["ClusterID"]))

def main(fasta_file, cluster_csv, gene_tsv, output_tsv):
    gene_to_cluster = load_clusters(cluster_csv)
    
    sequences_with_clusters = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        sequence = str(record.seq)
        cluster_id = gene_to_cluster.get(gene_id, "NA")
        sequences_with_clusters.append([gene_id, sequence, cluster_id])
    
    seq_clusters_df = pd.DataFrame(sequences_with_clusters, columns=["GeneID", "Sequence", "ClusterID"])
    
    gene_df = pd.read_csv(gene_tsv, sep="\t")
    gene_column = detect_sequence_column(gene_df, seq_clusters_df["Sequence"].values)
    
    sequence_to_cluster = dict(zip(seq_clusters_df["Sequence"], seq_clusters_df["ClusterID"]))
    gene_df["ClusterID"] = gene_df[gene_column].map(sequence_to_cluster).fillna("NA")
    gene_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Results saved to {output_tsv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map cluster IDs to sequences and merge with gene data.")
    parser.add_argument("--fasta_file", required=True, help="Input FASTA file with gene sequences")
    parser.add_argument("--cluster_csv", required=True, help="CSV file containing cluster information")
    parser.add_argument("--gene_tsv", required=True, help="TSV file containing gene data")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file path")

    args = parser.parse_args()
    main(args.fasta_file, args.cluster_csv, args.gene_tsv, args.output_tsv)