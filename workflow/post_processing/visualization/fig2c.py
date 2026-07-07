import pandas as pd
import matplotlib.pyplot as plt

# 1. Load data from a text file
# Replace 'data.txt' with your actual file path (e.g., 'input_data.tsv')
df = pd.read_csv('/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/04_filtered/sex_linked_species_lists.tsv', sep='\t')

# 2. Explode data and remove missing values
df['Sex_Linked_Genes'] = df['Sex_Linked_Genes'].fillna('')
df_genes = df.assign(Gene=df['Sex_Linked_Genes'].str.split(',')).explode('Gene')
df_genes = df_genes[df_genes['Gene'] != '']

# 3. Calculate gene frequency and sort (Core logic to group unique genes by Order)
gene_counts = df_genes.groupby('Gene')['Sample'].nunique()
gene_orders = df_genes.groupby('Gene')['Order'].first()

sort_df = pd.DataFrame({
    'Frequency': gene_counts,
    'Order': gene_orders
})

sort_df = sort_df.sort_values(by=['Frequency', 'Order'], ascending=[False, True])
ordered_genes = sort_df.index.tolist()

# 4. Y-axis (Order) mapping and apply user-provided colors
order_gene_df = df_genes[['Order', 'Gene']].drop_duplicates()
ordered_orders = sorted(order_gene_df['Order'].unique())
order_y_map = {order: i for i, order in enumerate(ordered_orders)}
order_gene_df['Y'] = order_gene_df['Order'].map(order_y_map)

# List of 12 requested colors (applied sequentially from the top)
user_colors = [
    '#E58700', '#C99800', '#A3A500', '#6BB100', 
    '#00BF7D', '#00C0AF', '#00BCD8', '#00B0F6', 
    '#619CFF', '#B983FF', '#E76BF3', '#FD61D1'
]

# To prevent errors if there are more than 12 Orders, colors are recycled from the beginning.
custom_colors = {order: user_colors[i % len(user_colors)] for i, order in enumerate(ordered_orders)}

# 5. Set up the plot canvas
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), gridspec_kw={'height_ratios': [1.2, 1]}, sharex=True)
plt.subplots_adjust(hspace=0.03)

# ---------------------------------------------------------
# [Top Plot] Bar Chart: Total samples per gene
# ---------------------------------------------------------
ax1.bar(ordered_genes, gene_counts[ordered_genes], color='#e0e0e0', edgecolor='#a0a0a0', width=0.7)
ax1.set_ylabel('Total Samples', fontsize=14, labelpad=15)
ax1.grid(axis='y', linestyle='--', alpha=0.5, color='lightgray')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

# ---------------------------------------------------------
# [Bottom Plot] Dot Plot: Gene distribution by Order
# ---------------------------------------------------------
# 1) Draw vertical lines
for gene in ordered_genes:
    gene_data = order_gene_df[order_gene_df['Gene'] == gene]
    if len(gene_data) > 1:
        ymin = gene_data['Y'].min()
        ymax = gene_data['Y'].max()
        ax2.vlines(x=gene, ymin=ymin, ymax=ymax, color='#d3d3d3', zorder=1, linewidth=1.2, linestyle=':')

# 2) Plot colored dots (mapped to the provided colors)
for order in ordered_orders:
    order_data = order_gene_df[order_gene_df['Order'] == order]
    c = custom_colors.get(order, '#95a5a6')
    ax2.scatter(order_data['Gene'], order_data['Y'], color=c, s=80, zorder=2, edgecolors='white', linewidths=0.5)

# Set Y-axis (Order) labels and invert the axis
ax2.set_yticks(range(len(ordered_orders)))
ax2.set_yticklabels(ordered_orders, fontsize=11, color='#333333')
ax2.invert_yaxis() 

# Rotate X-axis (Gene) labels by 90 degrees
ax2.set_xticks(range(len(ordered_genes)))
ax2.set_xticklabels(ordered_genes, rotation=90, fontsize=8, color='#555555')

# Clean up the bottom background
ax2.grid(axis='y', linestyle='-', alpha=0.4, color='lightgray')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False) 
ax2.spines['bottom'].set_visible(False)
ax2.tick_params(axis='y', length=0) 

# 6. Render and save the SVG file
plt.tight_layout()
output_filename = "/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/04_filtered/sex_linked_loci_distribution.svg"
plt.savefig(output_filename, format="svg", bbox_inches="tight", transparent=False, facecolor='white')
print(f"Figure successfully saved to: {output_filename}")