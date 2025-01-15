import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_performance(input_file, output_file):
    # Set global font size (all 20)
    plt.rcParams.update({'font.size': 20})
    
    # Read data
    df = pd.read_csv(input_file, sep='\t')
    
    # Use the Dataset Size column directly
    sizes = df['Dataset Size'].values
    scores = df['Manakov test set'].values
    
    # Create figure with white background
    plt.figure(figsize=(7, 6), facecolor='white')
    ax = plt.gca()
    ax.set_facecolor('white')
    
    # Plot data with exact styling
    plt.plot(sizes, scores, color='blue', linestyle='--', 
             marker='o', markersize=8, markerfacecolor='blue',
             markeredgecolor='blue', linewidth=1)
    
    # Configure axes
    plt.xscale('log')
    plt.xlim(10**2, 10**6)
    plt.ylim(0.4, 0.9)
    
    # Set grid style
    plt.grid(True, color='#E0E0E0', linestyle='-', alpha=0.5)
    
    # Configure ticks and labels with larger font
    plt.xticks([10**2, 10**3, 10**4, 10**5, 10**6], 
               ['2', '3', '4', '5', '6'],
               fontsize=20)
    plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
               ['0.40', '0.50', '0.60', '0.70', '0.80', '0.90'],
               fontsize=20)
    
    # Set labels with larger font
    plt.xlabel('log10(number of samples)', fontsize=20, labelpad=10)
    
    # Adjust spines
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')

def main():
    parser = argparse.ArgumentParser(description='Plot performance data')
    parser.add_argument('input', help='Input TSV file')
    parser.add_argument('output', help='Output PNG file')
    args = parser.parse_args()
    
    plot_performance(args.input, args.output)

if __name__ == '__main__':
    main()