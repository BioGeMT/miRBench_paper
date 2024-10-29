import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from seed_utils import find_seed_match

def process_dataset(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['seed_type'] = df.apply(lambda row: find_seed_match(row['gene'], row['noncodingRNA']), axis=1)
    
    grouped = df.groupby('noncodingRNA_fam')
    result = grouped.agg({
        'seed_type': lambda x: (x.isin(['Seed6mer', 'Seed7mer', 'Seed8mer']).sum() / len(x)) * 100,
        'noncodingRNA_fam': 'count'
    })
    
    result.columns = ['TotalCanonicalSeed%', 'EntryCount']
    result = result.reset_index().sort_values('TotalCanonicalSeed%')
    
    # Define the tick ranges
    tick_ranges = [(i, i+5) for i in range(0, 100, 5)]
    
    # Function to assign tick based on TotalCanonicalSeed%
    def assign_tick(percentage):
        for i, (lower, upper) in enumerate(tick_ranges):
            if lower <= percentage < upper:
                return i
        return 19  # For 100% (falls into the last tick)
    
    # Add the new 'ticks' column
    result['ticks'] = result['TotalCanonicalSeed%'].apply(assign_tick)
    
    return result

def plot_violin(dataframes, output_file):
    # Extract the file extension to determine the output format
    file_extension = os.path.splitext(output_file)[1].lower()
    
    fig, ax = plt.subplots(figsize=(3500 / 300, 3500 / 300), dpi=300)
    colors = ['#f0e442', '#0072b2', '#cc79a7']
    positions = np.arange(len(dataframes))

    for i, (dataset_name, df) in enumerate(dataframes.items()):
        # Count the number of miRNA families in each tick
        tick_counts = df['ticks'].value_counts().sort_index()
        
        # Create the violin plot data
        violin_data = []
        for tick, count in tick_counts.items():
            violin_data.extend([tick * 5 + 2.5] * count)  # Use middle of each 5% range
        
        parts = ax.violinplot(violin_data, positions=[i], showmeans=False, showmedians=False, showextrema=False, widths=0.7)  
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i % len(colors)])
            pc.set_edgecolor('black')
            pc.set_linewidth(1.5)  
            pc.set_alpha(1)

        median = np.median(violin_data)
        ax.hlines(median, i - 0.35, i + 0.35, color='black', linewidth=2.5)

    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 10))
    ax.set_yticklabels([f'{i}%' for i in range(0, 101, 10)], fontsize=28)
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_linewidth(1.5)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    ax.set_xticks(positions)
    ax.set_xticklabels(list(dataframes.keys()), fontsize=24, rotation=0, ha='center', weight='semibold')  
    ax.xaxis.set_tick_params(pad=15)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xlim(-0.5, len(positions) - 0.5)

    plt.tight_layout()
    
    # Save the plot based on the file extension
    if file_extension == '.svg':
        plt.savefig(output_file, format='svg', bbox_inches='tight')
    else:  # Default to PNG
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png')
    
    print(f"Plot saved as {output_file}")
    plt.close()

def create_output_tsv(dataframes, output_file):
    # Create a list to store data from all datasets
    all_data = []

    for dataset_name, df in dataframes.items():
        # Group by ticks and calculate the count of miRNA families in each tick
        tick_data = df.groupby('ticks').size().reset_index(name='Count')
        tick_data['Dataset'] = dataset_name
        tick_data['Percentage'] = tick_data['ticks'] * 5 + 2.5  # Middle of each 5% range
        all_data.append(tick_data)

    # Combine data from all datasets
    combined_data = pd.concat(all_data, ignore_index=True)

    # Sort the data first by Dataset, then by Percentage in descending order
    combined_data = combined_data.sort_values(['Dataset', 'Percentage'], ascending=[True, False])

    # Reorder columns and rename for clarity
    combined_data = combined_data[['Dataset', 'Percentage', 'Count']]
    combined_data.columns = ['Dataset', 'Canonical_Seed_Percentage', 'miRNA_Family_Count']

    # Save to TSV file
    combined_data.to_csv(output_file, sep='\t', index=False)
    print(f"Data saved to {output_file}")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process miRNA datasets, create violin plots, and output data to TSV.')
    parser.add_argument('input_files', nargs='+', help='Input TSV files')
    parser.add_argument('-o', '--output_plot', default='seed_prevalence_plot.png', 
                      help='Output file for the plot (supported formats: .png, .svg)')
    parser.add_argument('-t', '--output_tsv', default='seed_prevalence_data.tsv', 
                      help='Output TSV file for plot data')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    dataframes = {}
    for i, file_path in enumerate(args.input_files):
        dataset_name = os.path.splitext(os.path.basename(file_path))[0]
        df = process_dataset(file_path)
        
        # Filter the first dataset to include only miRNA families with 100 or more entries
        if i == 0:
            df = df[df['EntryCount'] >= 100]
        # Filter the second dataset to include only miRNA families with 10 or more entries
        elif i == 1:
            df = df[df['EntryCount'] >= 10]
        # Filter the third dataset to exclude miRNA families with 1 or fewer entries
        elif i == 2:
            df = df[df['EntryCount'] > 1]
        
        dataframes[dataset_name] = df
    
    plot_violin(dataframes, args.output_plot)
    create_output_tsv(dataframes, args.output_tsv)

if __name__ == "__main__":
    main()
