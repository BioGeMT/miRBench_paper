import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import argparse
from matplotlib.patches import Patch, Rectangle
import numpy as np

def read_data(file_path):
    return pd.read_csv(file_path, sep='\t')

def create_sets(df):
    manakov_set = set(df[df['Manakov2022'] > 0]['miRNA Family'])
    hejret_set = set(df[df['Hejret2023'] > 0]['miRNA Family'])
    klimentova_set = set(df[df['Klimentova2022'] > 0]['miRNA Family'])
    return manakov_set, hejret_set, klimentova_set

def create_venn_diagram(manakov_set, hejret_set, klimentova_set):
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Calculate all possible intersections
    subsets = (
        len(manakov_set - (hejret_set | klimentova_set)),                    # 100
        len(hejret_set - (manakov_set | klimentova_set)),                    # 010
        len(manakov_set & hejret_set - klimentova_set),                      # 110
        len(klimentova_set - (manakov_set | hejret_set)),                    # 001
        len(manakov_set & klimentova_set - hejret_set),                      # 101
        len(hejret_set & klimentova_set - manakov_set),                      # 011
        len(manakov_set & hejret_set & klimentova_set)                       # 111
    )
    
    v = venn3(subsets=subsets, set_labels=('Manakov2022', 'Hejret2023', 'Klimentova2022'))

    # Updated colors
    colors = ['#f0e442', '#0072b2', '#cc79a7']  # yellow, blue, purple

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
    
    # Original rectangle positions
    rectangle_positions = [
        (-0.68, 0.35),  # Manakov2022
        (0.55, 0.18),   # Hejret2023
        (0.52, -0.45)   # Klimentova2022
    ]

    # Add original rectangles
    for i, pos in enumerate(rectangle_positions):
        rect = Rectangle(pos, 0.05, 0.05, facecolor=colors[i], edgecolor='none')
        ax.add_artist(rect)

    # Set font size and weight for set labels
    for i, text in enumerate(v.set_labels):
        if text:
            text.set_fontsize(20)
            text.set_fontweight('medium')
            text.set_color('black')
            text.set_position(label_positions[i])

    # Dictionary to store original positions of labels
    original_positions = {}
    subset_texts = {}  # Store all subset label texts and positions
    
    # Set font size and weight for subset labels (numbers) and store positions
    for i, text in enumerate(v.subset_labels):
        if text:
            text.set_fontsize(20)
            text.set_fontweight('bold')
            text.set_color('black')
            # Store all visible numbers and their positions
            if text.get_text() != '0':  # Only store non-zero values
                subset_texts[f"{text.get_text()}_{i}"] = text.get_position()
            if text.get_text() in ['1', '5']:
                original_positions[text.get_text() + f"_{i}"] = text.get_position()
                text.set_visible(False)
            elif text.get_text() == '0':
                text.set_visible(False)

    # Add relocated labels and connecting lines
    relocated_x_position = 0.55
    
    for label_id, orig_pos in original_positions.items():
        label_text, i = label_id.split('_')
        i = int(i)
        new_pos = (relocated_x_position, orig_pos[1])
        
        # Draw straight horizontal connecting line
        ax.annotate('',
                   xy=orig_pos,
                   xytext=new_pos,
                   arrowprops=dict(arrowstyle='-',
                                 connectionstyle='arc3,rad=0',  # Straight line
                                 color='black',
                                 alpha=1.0,
                                 linewidth=2.5))
        
        # Add relocated label
        ax.text(new_pos[0], new_pos[1],
                label_text,
                fontsize=20,
                fontweight='medium',
                ha='left',
                va='center')
        
        rect_x = new_pos[0] + 0.03  # Position rectangles after the numbers
        
        # Add new rectangles based on specific regions
        
        # Blue rectangle next to Hejret unique (5)
        if label_text == '5' and i == 1:  # Hejret unique region
            rect_blue = Rectangle((rect_x, new_pos[1] - 0.025), 0.05, 0.05, 
                                facecolor=colors[1], edgecolor='none')
            ax.add_artist(rect_blue)
        
        # Yellow/Purple rectangle next to Manakov/Klimentova intersection (5)
        elif label_text == '5' and i == 4:  # Manakov/Klimentova intersection
            rect_yellow_half = Rectangle((rect_x, new_pos[1] - 0.025), 0.025, 0.05, 
                                      facecolor=colors[0], edgecolor='none')
            rect_purple_half = Rectangle((rect_x + 0.025, new_pos[1] - 0.025), 0.025, 0.05, 
                                      facecolor=colors[2], edgecolor='none')
            ax.add_artist(rect_yellow_half)
            ax.add_artist(rect_purple_half)
        
        # Blue/Purple rectangle next to Hejret/Klimentova intersection (1)
        elif label_text == '1' and i == 5:  # Hejret/Klimentova intersection
            rect_blue_half = Rectangle((rect_x, new_pos[1] - 0.025), 0.025, 0.05, 
                                    facecolor=colors[1], edgecolor='none')
            rect_purple_half2 = Rectangle((rect_x + 0.025, new_pos[1] - 0.025), 0.025, 0.05, 
                                       facecolor=colors[2], edgecolor='none')
            ax.add_artist(rect_blue_half)
            ax.add_artist(rect_purple_half2)

    # Remove axes
    plt.axis('off')
    return v

def save_diagram(output_path):
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def calculate_statistics(manakov_set, hejret_set, klimentova_set):
    # Calculate all possible intersections and their contents
    m_only = manakov_set - (hejret_set | klimentova_set)
    h_only = hejret_set - (manakov_set | klimentova_set)
    k_only = klimentova_set - (manakov_set | hejret_set)
    mh_only = manakov_set & hejret_set - klimentova_set
    mk_only = manakov_set & klimentova_set - hejret_set
    hk_only = hejret_set & klimentova_set - manakov_set
    mhk = manakov_set & hejret_set & klimentova_set
    
    # Basic statistics
    stats_data = {
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
            len(manakov_set | hejret_set | klimentova_set),
            len(mhk),
            len(m_only),
            len(h_only),
            len(k_only),
            len(manakov_set),
            len(hejret_set),
            len(klimentova_set)
        ]
    }
    
    # Detailed intersection statistics
    intersection_stats = {
        'Intersection': [
            'Manakov2022 only',
            'Hejret2023 only',
            'Klimentova2022 only',
            'Manakov2022 ∩ Hejret2023 only',
            'Manakov2022 ∩ Klimentova2022 only',
            'Hejret2023 ∩ Klimentova2022 only',
            'All three studies'
        ],
        'Count': [
            len(m_only),
            len(h_only),
            len(k_only),
            len(mh_only),
            len(mk_only),
            len(hk_only),
            len(mhk)
        ],
        'miRNA_Families': [
            ', '.join(sorted(m_only)) if len(m_only) < 20 else f"{len(m_only)} families",
            ', '.join(sorted(h_only)) if len(h_only) < 20 else f"{len(h_only)} families",
            ', '.join(sorted(k_only)) if len(k_only) < 20 else f"{len(k_only)} families",
            ', '.join(sorted(mh_only)) if len(mh_only) < 20 else f"{len(mh_only)} families",
            ', '.join(sorted(mk_only)) if len(mk_only) < 20 else f"{len(mk_only)} families",
            ', '.join(sorted(hk_only)) if len(hk_only) < 20 else f"{len(hk_only)} families",
            ', '.join(sorted(mhk)) if len(mhk) < 20 else f"{len(mhk)} families"
        ]
    }
    
    return pd.DataFrame(stats_data), pd.DataFrame(intersection_stats)

def main():
    parser = argparse.ArgumentParser(description='Create a Venn diagram from miRNA family data.')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file path')
    parser.add_argument('-o', '--output', required=True, help='Output image file path')
    parser.add_argument('-s', '--stats', required=True, help='Output statistics TSV file path')
    parser.add_argument('-d', '--detailed', required=True, help='Output detailed intersection statistics TSV file path')
    args = parser.parse_args()

    df = read_data(args.input)
    sets = create_sets(df)
    create_venn_diagram(*sets)
    save_diagram(args.output)
    
    stats_df, detailed_stats_df = calculate_statistics(*sets)
    stats_df.to_csv(args.stats, sep='\t', index=False)
    detailed_stats_df.to_csv(args.detailed, sep='\t', index=False)

if __name__ == "__main__":
    main()
