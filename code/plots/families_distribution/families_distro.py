import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import math
import numpy as np

def plot_mirna_families(input_file):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Sort the dataframe by Occurrences in descending order
    df_sorted = df.sort_values('Occurrences', ascending=False)

    # Get the top 100 families
    top_100 = df_sorted.head(100)

    # Calculate the sum of the rest
    rest_sum = df_sorted.iloc[100:]['Occurrences'].sum()

    # Create a new dataframe with top 100 and 'Other'
    other_df = pd.DataFrame({'miRNA Family': ['Other'], 'Occurrences': [rest_sum]})
    plot_data = pd.concat([top_100, other_df], ignore_index=True)

    # Set the style
    sns.set_style("white")

    # Determine the maximum y-axis value (round up to the nearest round value)
    max_value = plot_data['Occurrences'].max()
    y_max = 10 ** math.ceil(math.log10(max_value))

    # Get the first word of the input file name
    first_word = os.path.basename(input_file).split('_')[0]

    # Create the bar plot
    plt.figure(figsize=(80, 24))  # Increased figure size for better spacing

    # Adjust values to ensure visibility of columns with value 1
    adjusted_values = np.maximum(plot_data['Occurrences'], 1.1)  # Set minimum value to 1.1 for visibility
    ax = plt.bar(plot_data['miRNA Family'], adjusted_values, color='darkblue', width=0.7)

    # Customize the plot
    plt.title(f'Top 100 miRNA Families by Chimeric Reads - {first_word}', fontsize=40, fontweight='bold')
    plt.xlabel('miRNA Family', fontsize=32, fontweight='bold')
    plt.ylabel('Number of Chimeric Reads', fontsize=32, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=28, fontweight='bold')
    plt.yticks(fontsize=28, fontweight='bold')

    # Set y-axis to log scale
    plt.yscale('log')

    # Set custom y-axis ticks
    y_ticks = [1, 10, 100, 1000, 10000, 100000, 1000000]
    y_ticks = [tick for tick in y_ticks if tick <= y_max]  # Only include ticks up to y_max
    plt.yticks(y_ticks, [str(tick) for tick in y_ticks], fontsize=28, fontweight='bold')

    # Set y-axis limit
    plt.ylim(0.9, y_max)  # Start slightly below 1 to show full bar for value 1

    # Remove top and right spines
    sns.despine()

    # Adjust layout and increase bottom margin for x-axis labels
    plt.tight_layout(pad=4)
    plt.subplots_adjust(bottom=0.2)  # Increase bottom margin

    # Generate output file name
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}_families_distro.png"

    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')

    print(f"Plot has been saved as '{output_file}'")

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path_to_input_tsv_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    try:
        plot_mirna_families(input_file)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()