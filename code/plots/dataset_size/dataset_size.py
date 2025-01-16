import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_performance(input_file, output_file):

    plt.rcParams.update({'font.size': 20})
    
    # Read data
    df = pd.read_csv(input_file, sep='\t')
    
    # Filter out the 7720 dataset size
    df = df[df['Dataset Size'] != 7720]
 
    plt.figure(figsize=(7, 6), facecolor='white')
    ax = plt.gca()
    ax.set_facecolor('white')
    
    # Plot data with hollow markers
    plt.plot(df['Dataset Size'], df['Manakov test set'], color='blue', 
             marker='o', markersize=10, markerfacecolor='none',
             markeredgecolor='blue', linestyle='None', label='Manakov test set')
    plt.plot(df['Dataset Size'], df['Manakov leftout'], color='red', 
             marker='o', markersize=10, markerfacecolor='none',
             markeredgecolor='red', linestyle='None', label='Manakov leftout')
    plt.plot(df['Dataset Size'], df['Hejret new test set'], color='green', 
             marker='o', markersize=10, markerfacecolor='none',
             markeredgecolor='green', linestyle='None', label='Hejret new test set')
    plt.plot(df['Dataset Size'], df['Klimentova new test set'], color='purple', 
             marker='o', markersize=10, markerfacecolor='none',
             markeredgecolor='purple', linestyle='None', label='Klimentova new test set')
    
    # Configure axes
    plt.xscale('log')
    plt.xlim(10**2, 10**7)  # Extended to 10^7 to show full range
    plt.ylim(0.4, 0.9)
    
    # Set grid style
    plt.grid(True, color='#E0E0E0', linestyle='-', alpha=0.5)
    
    # Configure ticks and labels with larger font
    plt.xticks([10**2, 10**3, 10**4, 10**5, 10**6, 10**7], 
               ['2', '3', '4', '5', '6', '7'],
               fontsize=20)
    plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
               ['0.40', '0.50', '0.60', '0.70', '0.80', '0.90'],
               fontsize=20)
    
    # Set labels with larger font
    plt.xlabel('log10(dataset size)', fontsize=20, labelpad=10)
    
    # Add legend
    plt.legend(fontsize=16, loc='lower right')
    
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
