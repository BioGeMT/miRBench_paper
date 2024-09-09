import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import argparse

def read_data(file_path):
    return pd.read_csv(file_path, sep='\t')

def create_sets(df):
    manakov_set = set(df[df['Manakov2022'] > 0]['miRNA Family'])
    hejret_set = set(df[df['Hejret2023'] > 0]['miRNA Family'])
    klimentova_set = set(df[df['Klimentova2022'] > 0]['miRNA Family'])
    return manakov_set, hejret_set, klimentova_set

def create_venn_diagram(manakov_set, hejret_set, klimentova_set):
    plt.figure(figsize=(10, 10))
    v = venn3([manakov_set, hejret_set, klimentova_set], 
              set_labels=('Manakov2022', 'Hejret2023', 'Klimentova2022'))

    # set colors for each section
    colors = ['#f0e442', '#0072b2', '#009e73', '#e69f00', '#cc79a7', '#56b4e9', '#d55e00']
    for i, patch in enumerate(v.patches):
        if patch:
            patch.set_color(colors[i])
            patch.set_alpha(1)

    # draw circles
    venn3_circles([manakov_set, hejret_set, klimentova_set], linewidth=1.5, color='black')

    # set font size and weight for set labels
    for text in v.set_labels:
        text.set_fontsize(20)
        text.set_fontweight('bold')

    # set font size and weight for subset labels
    for text in v.subset_labels:
        if text:
            text.set_fontsize(16)
            text.set_fontweight('bold')

    # adjust position of set labels
    v.set_labels[0].set_position((-0.40, 0.50))
    v.set_labels[1].set_position((0.45, 0.25))
    v.set_labels[2].set_position((0.55, -0.35))

    # remove axes
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