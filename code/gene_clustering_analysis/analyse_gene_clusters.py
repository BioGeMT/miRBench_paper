import pandas as pd
import matplotlib.pyplot as plt
import argparse

def plot_cluster_size_distribution(df, plot_output):
    """
    Plots the distribution of cluster sizes.
    """
    # Count the occurrences for each cluster ID.
    cluster_counts = df["gene_cluster_ID"].value_counts()

    # Compute the distribution: number of clusters that have a given occurrence count.
    distribution = cluster_counts.value_counts().sort_index()

    plt.figure(figsize=(10, 6))
    distribution.plot(kind="bar")
    plt.xlabel("Cluster Size")
    plt.ylabel("Number of Clusters")
    plt.title("Distribution of Cluster Sizes")
    plt.tight_layout()
    plt.savefig(plot_output)
    plt.close()

def compute_biggest_cluster_percentage(df):
    """
    Computes the percentage of total occurrences that are in the biggest cluster.
    """
    # Count occurrences per cluster.
    cluster_counts = df["gene_cluster_ID"].value_counts()
    
    # Identify the cluster with the highest occurrence count.
    biggest_cluster_id = cluster_counts.idxmax()
    biggest_cluster_count = cluster_counts.max()
    
    # Total number of occurrences (rows) in the DataFrame.
    total_occurrences = df.shape[0]
    
    # Calculate the percentage of occurrences in the biggest cluster.
    percentage = (biggest_cluster_count / total_occurrences) * 100
    
    return biggest_cluster_id, biggest_cluster_count, percentage

def compute_singleton_cluster_percentage(df):
    """
    Computes the percentage of total occurrences that are in singleton clusters.
    """
    # Count occurrences per cluster.
    cluster_counts = df["gene_cluster_ID"].value_counts()
    
    # Identify singleton clusters: clusters with exactly one occurrence.
    singleton_clusters = cluster_counts[cluster_counts == 1]
    
    # Each singleton cluster contributes one occurrence.
    singleton_occurrences = singleton_clusters.sum()
    
    total_occurrences = df.shape[0]
    percentage = (singleton_occurrences / total_occurrences) * 100
    
    return percentage

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to the input file (tab-separated)")
    parser.add_argument("--plot_output", required=True, help="Path to the output plot file (.png)")
    parser.add_argument("--metrics_output", required=True, help="Path to the output metrics file (.txt)")
    return parser.parse_args()

def main():
    args = parse_args()
    # Read the input file (assumed to be tab-separated).
    df = pd.read_csv(args.input, sep="\t")
    
    # 1. Plot the distribution of cluster sizes.
    plot_cluster_size_distribution(df, args.plot_output)
    
    # 2. Compute metrics for the biggest cluster.
    biggest_cluster_id, biggest_cluster_count, biggest_cluster_pct = compute_biggest_cluster_percentage(df)
    
    # 3. Compute the percentage of occurrences in singleton clusters.
    singleton_cluster_pct = compute_singleton_cluster_percentage(df)
    
    # Write the computed metrics to an output file.
    with open(args.metrics_output, "w") as f:
        f.write(f"Biggest Cluster ID: {biggest_cluster_id}\n")
        f.write(f"Count of Occurrences in Biggest Cluster: {biggest_cluster_count}\n")
        f.write(f"% of Occurrences in Biggest Cluster: {biggest_cluster_pct:.2f}%\n")
        f.write(f"% of Occurrences in Singleton Clusters: {singleton_cluster_pct:.2f}%\n")

if __name__ == "__main__":
    main()
