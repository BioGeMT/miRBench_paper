import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

def load_data(file_path):
    try:
        return pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

def process_data(data):
    categories = ["5' UTR", "3' UTR", 'intron', 'exon', 'Unknown']
    
    def categorize_feature(feature):
        feature = str(feature).lower()
        if 'five_prime_utr' in feature:
            return "5' UTR"
        elif 'three_prime_utr' in feature:
            return "3' UTR"
        elif 'intron' in feature:
            return 'intron'
        elif 'exon' in feature:
            return 'exon'
        else:
            return 'Unknown'

    data['category'] = data['feature'].apply(categorize_feature)
    
    # Count total occurrences for each miRNA family
    total_occurrences = data['noncodingRNA_fam'].value_counts()
    
    # Group by miRNA family and category, then count occurrences
    grouped = data.groupby(['noncodingRNA_fam', 'category']).size().unstack(fill_value=0)
    
    for cat in categories:
        if cat not in grouped.columns:
            grouped[cat] = 0
    
    grouped = grouped[categories]
    
    # Calculate percentages
    total = grouped.sum(axis=1)
    percentages = grouped.div(total, axis=0) * 100
    
    # Add total occurrences to the percentages dataframe
    percentages['Total_Occurrences'] = total_occurrences
    
    # Sort by total occurrences and get top 100
    top_100 = percentages.nlargest(100, 'Total_Occurrences')
    
    # Remove the Total_Occurrences column for plotting
    top_100 = top_100.drop('Total_Occurrences', axis=1)
    
    return top_100

def create_plot(data, input_file):
    # Set the style
    sns.set_style("white")

    # Get the first word of the input file name
    first_word = input_file.split('_')[0]

    # Create the plot
    fig, ax = plt.subplots(figsize=(80, 24))

    # Define colors for each feature
    colors = ['#FFA07A', '#98FB98', '#87CEFA', '#DDA0DD', 'gray']

    # Create the stacked bar chart
    data.plot(kind='bar', stacked=True, ax=ax, color=colors)

    # Customize the plot
    plt.title(f'Top 100 miRNA Families by Chimeric Reads - {first_word}', fontsize=40, fontweight='bold')
    plt.xlabel('miRNA Family', fontsize=32, fontweight='bold')
    plt.ylabel('Percentage of Chimeric Reads', fontsize=32, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=28, fontweight='bold')
    plt.yticks(fontsize=28, fontweight='bold')

    # Set y-axis to percentage
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.set_yticklabels([f'{x}%' for x in range(0, 101, 20)], fontsize=28, fontweight='bold')

    # Add grid lines
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)

    # Remove top and right spines
    sns.despine()

    # Adjust layout and increase bottom margin for x-axis labels
    plt.tight_layout(pad=4)
    plt.subplots_adjust(bottom=0.2)  # Increase bottom margin

    # Adjust legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=24, frameon=False)

    return fig

def main():
    parser = argparse.ArgumentParser(description="miRNA Family Distribution Plot")
    parser.add_argument('--ifile', help="Input file (TSV format)", required=True)
    parser.add_argument('--ofile', help="Output file for plot (e.g., plot.png)", required=True)
    args = parser.parse_args()

    # Load data
    data = load_data(args.ifile)

    # Process data
    processed_data = process_data(data)

    # Create plot
    fig = create_plot(processed_data, args.ifile)

    # Save plot
    fig.savefig(args.ofile, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {args.ofile}")

if __name__ == "__main__":
    main()