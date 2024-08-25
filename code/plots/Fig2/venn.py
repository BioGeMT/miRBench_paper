import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

# Read the TSV file
df = pd.read_csv('all_mirna_families_counts_sorted.tsv', sep='\t')

# Create sets for each study
manakov_set = set(df[df['Manakov2022'] > 0]['miRNA Family'])
hejret_set = set(df[df['Hejret2023'] > 0]['miRNA Family'])
klimentova_set = set(df[df['Klimentova2022'] > 0]['miRNA Family'])

# Set up the figure
plt.figure(figsize=(10, 10))

# Create the Venn diagram with the specified colors
v = venn3([manakov_set, hejret_set, klimentova_set], 
          set_labels=('Manakov2022', 'Hejret2023', 'Klimentova2022'))

# Assign specific colors to each region
colors = [
    '#f0e442',  # Set 1: Manakov2022 (Bright Yellow)
    '#0072b2',  # Set 2: Hejret2023 (Deep Sky Blue)
    '#009e73',  # Set 3: Klimentova2022 (Light Purple)
    '#cc79a7',  # Set 1 & 2 (Teal Green)
    '#e69f00',  # Set 1 & 3 (Golden Yellow)
    '#56b4e9',  # Set 2 & 3 (Soft Blue)
    '#d55e00'   # Set all (Burnt Orange)
]
for i, patch in enumerate(v.patches):
    if patch:  # Check if patch is not None (some regions might not exist)
        patch.set_color(colors[i])
        patch.set_alpha(1)  # Set full opacity

# Add circle outlines
venn3_circles([manakov_set, hejret_set, klimentova_set], linewidth=1.5, color='black')

# Customize text
for text in v.set_labels:
    text.set_fontsize(20)
    text.set_fontweight('bold')

for text in v.subset_labels:
    if text:
        text.set_fontsize(16)
        text.set_fontweight('bold')
        

# Adjust the position of the dataset names (distance from circles)
v.set_labels[0].set_position((-0.40, 0.50))  # Manakov2022 label position
v.set_labels[1].set_position((0.63, 0.15))   # Hejret2023 label position
v.set_labels[2].set_position((0, -0.60))     # Klimentova2022 label position

# Remove axes
plt.axis('off')

# Save the figure
plt.tight_layout()
plt.savefig('mirna_venn_diagram_proportional1.png', dpi=300, bbox_inches='tight')
plt.close()

# Calculate and print statistics
total_unique = len(manakov_set | hejret_set | klimentova_set)
common_to_all = len(manakov_set & hejret_set & klimentova_set)

print(f"Total unique miRNA families: {total_unique}")
print(f"miRNA families common to all studies: {common_to_all}")
print(f"Manakov2022 unique: {len(manakov_set - (hejret_set | klimentova_set))}")
print(f"Hejret2023 unique: {len(hejret_set - (manakov_set | klimentova_set))}")
print(f"Klimentova2022 unique: {len(klimentova_set - (manakov_set | hejret_set))}")
print(f"Manakov2022 total: {len(manakov_set)}")
print(f"Hejret2023 total: {len(hejret_set)}")
print(f"Klimentova2022 total: {len(klimentova_set)}")
