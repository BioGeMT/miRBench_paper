import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_performance(input_file, output_file):
    plt.rcParams.update({'font.size': 16})  # Reduced from 20
    
    # Read data
    df = pd.read_csv(input_file, sep='\t')
    
    # Filter out the 7720 dataset size
    df = df[df['Dataset Size'] != 7720]
    
    plt.figure(figsize=(7, 6), facecolor='white')
    ax = plt.gca()
    ax.set_facecolor('white')
    
    # Plot data with solid lines and hollow markers
    plt.plot(df['Dataset Size'], df['Manakov test set'], color='#0072b2', 
             marker='o', markersize=8, markerfacecolor='white',
             markeredgecolor='#0072b2', markeredgewidth=2,
             linestyle='-', linewidth=1.5, label='Manakov2022 test')
    
    plt.plot(df['Dataset Size'], df['Manakov leftout'], color='#56b4e9', 
             marker='s', markersize=8, markerfacecolor='white',
             markeredgecolor='#56b4e9', markeredgewidth=2,
             linestyle='-', linewidth=1.5, label='Manakov2022 left-out')
    
    plt.plot(df['Dataset Size'], df['Hejret new test set'], color='#009e73', 
             marker='^', markersize=8, markerfacecolor='white',
             markeredgecolor='#009e73', markeredgewidth=2,
             linestyle='-', linewidth=1.5, label='Hejret2023 test')
    
    plt.plot(df['Dataset Size'], df['Klimentova new test set'], color='#d55e00', 
             marker='D', markersize=8, markerfacecolor='white',
             markeredgecolor='#d55e00', markeredgewidth=2,
             linestyle='-', linewidth=1.5, label='Klimentova2022 test')
    
    # Configure axes
    plt.xscale('log')
    plt.xlim(10**2, 10**7)
    plt.ylim(0.4, 0.9)
    
    # Set grid style
    plt.grid(True, color='#E0E0E0', linestyle='-', alpha=0.5)
    
    # Configure ticks and labels
    plt.xticks([10**2, 10**3, 10**4, 10**5, 10**6, 10**7], 
               ['2', '3', '4', '5', '6', '7'],
               fontsize=16)
    plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
               ['0.40', '0.50', '0.60', '0.70', '0.80', '0.90'],
               fontsize=16)
    
    # Set labels
    plt.xlabel('log10(dataset size)', fontsize=16, labelpad=10)
    plt.ylabel('APS', fontsize=16, labelpad=10)
    
    # Add legend with smaller font
    plt.legend(fontsize=12, loc='lower right', 
              handlelength=3, handletextpad=0.5, 
              borderpad=0.4, borderaxespad=0.5)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust remaining spines
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(1.0)
    
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