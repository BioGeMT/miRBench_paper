import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data
df1 = pd.read_csv('top_10_mirna_families_dataset1.tsv', sep='\t')
df2 = pd.read_csv('top_10_mirna_families_dataset2.tsv', sep='\t')
df3 = pd.read_csv('top_10_mirna_families_dataset3.tsv', sep='\t')

# Specific names for the legend
legend_names = ['Manakov2022', 'Hejret2023', 'Klimentova2022']

# Function to calculate percentages
def calculate_percentages(df):
    total = df['Count'].sum()
    df['Percentage'] = df['Count'] / total * 100
    return df

# Calculate percentages for each dataset
df1 = calculate_percentages(df1)
df2 = calculate_percentages(df2)
df3 = calculate_percentages(df3)

# Combine all datasets
all_data = pd.concat([df1, df2, df3])

# Group by miRNA Family and sum the counts
grouped_data = all_data.groupby('miRNA Family')['Count'].sum().reset_index()

# Get top 10 overall families
top_10_overall = grouped_data.sort_values('Count', ascending=False).head(10)['miRNA Family'].tolist()

# Combine all unique families while preserving the order
all_families = top_10_overall 

# Function to get percentage for a family (or 0 if not present)
def get_percentage(df, family):
    if family in df['miRNA Family'].values:
        return df[df['miRNA Family'] == family]['Percentage'].values[0]
    return 0

# Prepare data for plotting
x = np.arange(len(all_families))
width = 0.25

percentages1 = [get_percentage(df1, family) for family in all_families]
percentages2 = [get_percentage(df2, family) for family in all_families]
percentages3 = [get_percentage(df3, family) for family in all_families]

# Set up the color palette (first 3 colors as specified)
colors = ['#f0e442', '#0072b2', '#cc79a7']

# Create the plot with adjusted figure size
fig, ax = plt.subplots(figsize=(30, 10), dpi=300)  # Increased width to accommodate more x-axis labels

# Function to plot bars centered around the tick
def plot_centered_bars(ax, x, percentages1, percentages2, percentages3, width, colors, legend_names):
    rects = []
    x_centers = []  # To store the center position of each group
    for i, (p1, p2, p3) in enumerate(zip(percentages1, percentages2, percentages3)):
        non_zero = [p for p in [p1, p2, p3] if p > 0]
        n_bars = len(non_zero)
        
        if n_bars == 0:
            x_centers.append(x[i])
            continue
        
        total_width = n_bars * width
        start = x[i] - total_width / 2 + width / 2
        
        center = start + (total_width - width) / 2  # Calculate center of the group
        x_centers.append(center)
        
        if p1 > 0:
            rect = ax.bar(start, p1, width, label=legend_names[0] if i == 0 else "", color=colors[0], edgecolor='black', linewidth=0.5)
            rects.extend(rect)
            start += width
        if p2 > 0:
            rect = ax.bar(start, p2, width, label=legend_names[1] if i == 0 else "", color=colors[1], edgecolor='black', linewidth=0.5)
            rects.extend(rect)
            start += width
        if p3 > 0:
            rect = ax.bar(start, p3, width, label=legend_names[2] if i == 0 else "", color=colors[2], edgecolor='black', linewidth=0.5)
            rects.extend(rect)
    
    return rects, x_centers

rects, x_centers = plot_centered_bars(ax, x, percentages1, percentages2, percentages3, width, colors, legend_names)

# Remove labels and title
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_title('')

# Set x-ticks to the center of each group
ax.set_xticks(x_centers)
ax.set_xticklabels(all_families, rotation=0, ha='center', fontsize=34, weight='semibold')

# Adjust y-axis to 50% max
ax.set_ylim(0, 50)
ax.set_yticks(range(0, 51, 10))
ax.set_yticklabels([f'{i}%' for i in range(0, 51, 10)], fontsize=34)
ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

# Increase width of x and y axis lines
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Move x-axis ticks down
ax.xaxis.set_tick_params(pad=10)

# Add legend without border and with semibold text
legend = ax.legend(fontsize=34, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)

# Make legend text semibold
for text in legend.get_texts():
    text.set_fontweight('semibold')

plt.tight_layout()

# Save the plot as a PNG file
plt.savefig('output.png', dpi=300, bbox_inches='tight')
print("Plot saved as output.png")