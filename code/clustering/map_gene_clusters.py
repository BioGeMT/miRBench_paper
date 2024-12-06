import argparse
import pandas as pd

def parse_fasta(fasta_path):

    seq_to_id = {}
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
            else:
                seq_to_id[line] = current_id
    return seq_to_id

def load_clusters(cluster_csv):
    
    clusters_df = pd.read_csv(cluster_csv)
    return dict(zip(clusters_df["Seq_ID"], clusters_df["Cluster_ID"]))

def main(cluster_csv, dataset_tsv, output_tsv, fasta_file):
  
    gene_to_cluster = load_clusters(cluster_csv)
    seq_to_id = parse_fasta(fasta_file)

    gene_df = pd.read_csv(dataset_tsv, sep="\t")
    gene_df["gene_cluster_ID"] = gene_df["gene"].map(
        lambda seq: gene_to_cluster.get(seq_to_id.get(seq))
    )

    gene_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Results saved to {output_tsv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map cluster IDs to sequences and merge with gene data.")
    parser.add_argument("--cluster_csv", required=True, help="CSV file containing cluster information")
    parser.add_argument("--dataset_tsv", required=True, help="Input dataset TSV file")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file path")
    parser.add_argument("--genes_fasta", required=True, help="Input FASTA file with Seq_IDs")

    args = parser.parse_args()
    main(args.cluster_csv, args.dataset_tsv, args.output_tsv, args.genes_fasta)

