import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
from collections import Counter
import os

def load_data(file_path):
    try:
        if file_path:
            return pd.read_csv(file_path, sep='\t')
        else:
            return pd.read_csv(sys.stdin, sep='\t')
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

def categorize_features(data):
    categories = ['five_prime_utr', 'three_prime_utr', 'intron', 'exon', 'Unknown']
    counts = {cat: 0 for cat in categories}
    other_features = Counter()
    total_entries = len(data)
    
    for feature in data['feature']:
        feature = str(feature).lower()
        if 'five_prime_utr' in feature:
            counts['five_prime_utr'] += 1
        elif 'three_prime_utr' in feature:
            counts['three_prime_utr'] += 1
        elif 'intron' in feature:
            counts['intron'] += 1
        elif 'exon' in feature:
            counts['exon'] += 1
        else:
            counts['Unknown'] += 1
            other_features[feature] += 1
    
    percentages = {cat: (count / total_entries) * 100 for cat, count in counts.items()}
    return percentages, other_features

def create_plot(percentages1, percentages2, percentages3, label1, label2, label3):
    categories = ["5' UTR", "3' UTR", 'Intron', 'Exon', 'Unknown']
    dataset1_values = [percentages1['five_prime_utr'], percentages1['three_prime_utr'], percentages1['intron'], percentages1['exon'], percentages1['Unknown']]
    dataset2_values = [percentages2['five_prime_utr'], percentages2['three_prime_utr'], percentages2['intron'], percentages2['exon'], percentages2['Unknown']]
    dataset3_values = [percentages3['five_prime_utr'], percentages3['three_prime_utr'], percentages3['intron'], percentages3['exon'], percentages3['Unknown']]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    bottom1, bottom2, bottom3 = 0, 0, 0
    colors = ['#FFA07A', '#98FB98', '#87CEFA', '#DDA0DD', 'gray']
    
    bars = []
    for i, cat in enumerate(categories):
        bar1 = ax.bar(label1, dataset1_values[i], bottom=bottom1, color=colors[i])
        ax.bar(label2, dataset2_values[i], bottom=bottom2, color=colors[i])
        ax.bar(label3, dataset3_values[i], bottom=bottom3, color=colors[i])
        bottom1 += dataset1_values[i]
        bottom2 += dataset2_values[i]
        bottom3 += dataset3_values[i]
        bars.append(bar1[0])
    
    ax.set_ylabel('Percentage')
    ax.set_title('miRNA Feature Distribution')
    ax.legend(bars, categories, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    return fig

def print_other_features(other_features, dataset_name):
    print(f"\nDetailed breakdown of 'Unknown' features for {dataset_name}:")
    total = sum(other_features.values())
    sorted_features = sorted(other_features.items(), key=lambda x: x[1], reverse=True)
    for feature, count in sorted_features:
        percentage = (count / total) * 100
        print(f"  {feature}: {count} ({percentage:.2f}%)")

def main():
    parser = argparse.ArgumentParser(description="miRNA Feature Distribution Plot")
    parser.add_argument('--ifile1', help="Input file for Dataset 1 (default: STDIN)")
    parser.add_argument('--ifile2', help="Input file for Dataset 2")
    parser.add_argument('--ifile3', help="Input file for Dataset 3")
    parser.add_argument('--ofile', help="Output file for plot (e.g., plot.png)")
    args = parser.parse_args()

    # Load datasets
    dataset1 = load_data(args.ifile1)
    dataset2 = load_data(args.ifile2)
    dataset3 = load_data(args.ifile3)

    # Get file names for labels
    label1 = os.path.basename(args.ifile1) if args.ifile1 else "STDIN"
    label2 = os.path.basename(args.ifile2) if args.ifile2 else "Dataset2"
    label3 = os.path.basename(args.ifile3) if args.ifile3 else "Dataset3"

    # Calculate percentages and get other features
    percentages1, other_features1 = categorize_features(dataset1)
    percentages2, other_features2 = categorize_features(dataset2)
    percentages3, other_features3 = categorize_features(dataset3)

    # Create plot
    fig = create_plot(percentages1, percentages2, percentages3, label1, label2, label3)

    # Save plot if output file is specified
    if args.ofile:
        fig.savefig(args.ofile, bbox_inches='tight')
    else:
        plt.show()

    # Print percentages and other features for verification
    print(f"{label1} percentages:", percentages1)
    print(f"{label2} percentages:", percentages2)
    print(f"{label3} percentages:", percentages3)
    
    print_other_features(other_features1, label1)
    print_other_features(other_features2, label2)
    print_other_features(other_features3, label3)

if __name__ == "__main__":
    main()