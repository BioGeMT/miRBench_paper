import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import argparse
from matplotlib.patches import Patch, Rectangle

def read_data(file_path):
    return pd.read_csv(file_path, sep='\t')

def create_sets(df):
    manakov_set = set(df[df['Manakov2022'] > 0]['miRNA Family'])
    hejret_set = set(df[df['Hejret2023'] > 0]['miRNA Family'])
    klimentova_set = set(df[df['Klimentova2022'] > 0]['miRNA Family'])
    return manakov_set, hejret_set, klimentova_set

def create_venn_diagram(manakov_set, hejret_set, klimentova_set):
    fig, ax = plt.subplots(figsize=(12, 10))
    v = venn3([manakov_set, hejret_set, klimentova_set], 
              set_labels=('Manakov2022', 'Hejret2023', 'Klimentova2022'))

    # Updated colors
    colors = ['#f0e442', '#0072b2', '#cc79a7']

    # Remove fill colors
    for patch in v.patches:
        if patch:
            patch.set_alpha(0)

    # Define line width for colored circles
    color_width = 10.0

    # Draw colored circles
    c_color = venn3_circles([manakov_set, hejret_set, klimentova_set], linewidth=color_width)
    
    for i, color_circle in enumerate(c_color):
        color_circle.set_edgecolor(colors[i])
        color_circle.set_zorder(1)

    # Configurable positions for labels and rectangles
    label_positions = [
        (-0.50, 0.42),  # Manakov2022
        (0.45, 0.25),   # Hejret2023
        (0.55, -0.35)   # Klimentova2022
    ]
    
    rectangle_positions = [
        (-0.68, 0.35),  # Manakov2022
        (0.55, 0.18),   # Hejret2023
        (0.52, -0.45)   # Klimentova2022
    ]

    # Set font size and weight for set labels
    for i, text in enumerate(v.set_labels):
        text.set_fontsize(20)
        text.set_fontweight('bold')
        text.set_color('black')
        text.set_position(label_positions[i])

    # Add rectangles
    for i, pos in enumerate(rectangle_positions):
        rect = Rectangle(pos, 0.05, 0.05, facecolor=colors[i], edgecolor='none')
        ax.add_artist(rect)

    # Set font size and weight for subset labels (numbers) and adjust positions
    for i, text in enumerate(v.subset_labels):
        if text:
            text.set_fontsize(15)
            text.set_fontweight('bold')
            text.set_color('black')
            if text.get_text() == '0':
                text.set_visible(False)
            else:
                # Adjust position of non-zero labels
                position = text.get_position()
                if i == 0:  # '100' label (Manakov2022 unique)
                    text.set_position((position[0] - 0.15, position[1] + 0.08))
                elif i == 1:  # '010' label (Hejret2023 unique)
                    text.set_position((position[0] - 0.00, position[1] + 0.0))
                elif i == 2:  # '110' label (Manakov2022 & Hejret2023)
                    text.set_position((position[0], position[1] + 0.08))
                elif i == 3:  # '001' label (Klimentova2022 unique)
                    text.set_position((position[0], position[1] - 0.08))
                elif i == 4:  # '101' label (Manakov2022 & Klimentova2022)
                    text.set_position((position[0] + 0.005, position[1] + 0.005))
                elif i == 5:  # '011' label (Hejret2023 & Klimentova2022)
                    text.set_position((position[0] + 0.05, position[1] - 0.02))
                elif i == 6:  # '111' label (common to all)
                    text.set_position((position[0], position[1] + 0.02))

    # Remove axes
    plt.axis('off')
    return v

def save_diagram(output_path):
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def calculate_statistics(manakov_set, hejret_set, klimentova_set):
    total_unique = len(manakov_set | hejret_set | klimentova_set)
    common_to_all = len(manakov_set & hejret_set & klimentova_set)
    manakov_unique = len(manakov_set - (hejret_set | klimentova_set))
    hejret_unique = len(hejret_set - (manakov_set | klimentova_set))
    klimentova_unique = len(klimentova_set - (manakov_set | hejret_set))
    manakov_total = len(manakov_set)
    hejret_total = len(hejret_set)
    klimentova_total = len(klimentova_set)

    return pd.DataFrame({
        'Statistic': [
            'Total unique miRNA families',
            'miRNA families common to all studies',
            'Manakov2022 unique',
            'Hejret2023 unique',
            'Klimentova2022 unique',
            'Manakov2022 total',
            'Hejret2023 total',
            'Klimentova2022 total'
        ],
        'Value': [
            total_unique,
            common_to_all,
            manakov_unique,
            hejret_unique,
            klimentova_unique,
            manakov_total,
            hejret_total,
            klimentova_total
        ]
    })

def save_statistics(stats_df, output_path):
    stats_df.to_csv(output_path, sep='\t', index=False)
    print(f"Statistics have been saved to {output_path}")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create a Venn diagram from miRNA family data.')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file path')
    parser.add_argument('-o', '--output', required=True, help='Output image file path')
    parser.add_argument('-s', '--stats', required=True, help='Output statistics TSV file path')
    return parser.parse_args()

def process_data(input_file):
    df = read_data(input_file)
    return create_sets(df)

def create_and_save_diagram(sets, output_file):
    v = create_venn_diagram(*sets)
    save_diagram(output_file)

def generate_and_save_statistics(sets, stats_file):
    stats_df = calculate_statistics(*sets)
    save_statistics(stats_df, stats_file)

def main():
    args = parse_arguments()
    sets = process_data(args.input)
    create_and_save_diagram(sets, args.output)
    generate_and_save_statistics(sets, args.stats)

if __name__ == "__main__":
    main()
