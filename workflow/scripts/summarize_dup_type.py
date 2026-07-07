import pandas as pd
import sys
import os
import argparse

# ---------------------------------------------------------
# 0. Argument Parsing
# ---------------------------------------------------------
parser = argparse.ArgumentParser(description="Summarize classified duplications by parsing Order from mapped_chromosomes.")
parser.add_argument("--classified_tsvs", nargs='+', required=True, help="List of classified TSV files")
parser.add_argument("--summary_by_order", required=True, help="Output path for summary by order")
parser.add_argument("--combined_all", required=True, help="Output path for combined classified results")

args = parser.parse_args()

input_tsvs = args.classified_tsvs
out_summary = args.summary_by_order
out_combined = args.combined_all

# ---------------------------------------------------------
# 1. Combine all classified results
# ---------------------------------------------------------
all_dfs = []

for tsv in input_tsvs:
    if os.stat(tsv).st_size == 0:
        continue
        
    df = pd.read_csv(tsv, sep='\t')
    if df.empty:
        continue
        
    # Extract sample ID from filename
    basename = os.path.basename(tsv)
    sample_id = basename.replace('_classified.tsv', '')
    df['sample_id'] = sample_id
    all_dfs.append(df)

if not all_dfs:
    # If no data at all, create empty output files
    open(out_summary, 'w').close()
    open(out_combined, 'w').close()
    sys.exit(0)

combined_df = pd.concat(all_dfs, ignore_index=True)

# ---------------------------------------------------------
# 2. Parse 'Order' and 'Gene_ID'
# ---------------------------------------------------------
# Example mapped_chromosomes: S0001_Polytrichales_Atrichum_angustatum_Chr1,S0001...
# 1) Order: Take the first hit (split by ','), split by '_', take index 1 -> 'Polytrichales'
# 2) Gene_ID: From qseqid (e.g., G4519|polytrichales...), split by '|', take index 0 -> 'G4519'

def parse_order(mapped_chromosomes):
    if pd.isna(mapped_chromosomes): return 'Unknown'
    # Multiple hits are separated by commas, so we take the first one
    first_chrom = str(mapped_chromosomes).split(',')[0]
    parts = first_chrom.split('_')
    
    # Extract the second element (index 1) which represents the Order
    if len(parts) > 1:
        return parts[1].capitalize() # Standardize format (e.g., 'Polytrichales')
    return 'Unknown'

def parse_gene_id(qseqid):
    if pd.isna(qseqid): return 'Unknown'
    return str(qseqid).split('|')[0]

# Apply parsing functions
combined_df['Order'] = combined_df['mapped_chromosomes'].apply(parse_order)
combined_df['Gene_ID'] = combined_df['qseqid'].apply(parse_gene_id)

# Export combined data
combined_df.to_csv(out_combined, sep='\t', index=False)

# ---------------------------------------------------------
# 3. Generate the Summary by Order
# ---------------------------------------------------------
summary_records = []

for order, group in combined_df.groupby('Order'):
    # Number of unique samples in this Order that have AT LEAST one duplication
    total_samples = group['sample_id'].nunique()
    
    # Count ALL types of duplications for this Order (Matching exactly with classify step)
    type_counts = group['duplication_type'].value_counts()
    
    # Autosomal / Unknown counts
    local_count = type_counts.get('Local_Duplication', 0)
    intra_count = type_counts.get('Intra_chromosomal_Duplication', 0)
    inter_chrom_count = type_counts.get('Inter_Chromosomal_Duplication', 0)
    pot_inter_chrom_count = type_counts.get('Potential_Inter_Chrom_Dups', 0)
    
    # Sex Chromosome specific counts
    local_sex_count = type_counts.get('Local_Sex_Chrom_Duplication', 0)
    intra_sex_count = type_counts.get('Intra_Sex_Chrom_Duplication', 0)
    inter_sex_count = type_counts.get('Inter_Sex_Chrom_Duplication', 0)
    sex_linked_count = type_counts.get('Sex_Linked_Duplication', 0)
    
    # Unclassified
    unclassified_count = type_counts.get('Unclassified_Complex', 0)

    # Identify highly conserved duplicated genes (>50% of samples in this Order)
    threshold = total_samples * 0.5
    
    # Conserved Inter-Chromosomal (WGD Candidate)
    wgd_events = group[group['duplication_type'] == 'Inter_Chromosomal_Duplication']
    if not wgd_events.empty:
        # Group by sample_id and Gene_ID to prevent overcounting if duplicated multiple times in one sample
        wgd_unique_per_sample = wgd_events.drop_duplicates(subset=['sample_id', 'Gene_ID'])
        wgd_gene_counts = wgd_unique_per_sample['Gene_ID'].value_counts()
        
        conserved_wgd = wgd_gene_counts[wgd_gene_counts >= threshold].index.tolist()
        conserved_wgd_str = ",".join(conserved_wgd) if conserved_wgd else "-"
    else:
        conserved_wgd_str = "-"
        
    # Conserved Sex-Linked
    sex_linked_events = group[group['duplication_type'] == 'Sex_Linked_Duplication']
    if not sex_linked_events.empty:
        sex_unique_per_sample = sex_linked_events.drop_duplicates(subset=['sample_id', 'Gene_ID'])
        sex_gene_counts = sex_unique_per_sample['Gene_ID'].value_counts()
        
        conserved_sex = sex_gene_counts[sex_gene_counts >= threshold].index.tolist()
        conserved_sex_str = ",".join(conserved_sex) if conserved_sex else "-"
    else:
        conserved_sex_str = "-"

    # Append to summary list
    summary_records.append({
        'Order': order,
        'Total_Samples_With_Dups': total_samples,
        'Local_Dups': local_count,
        'Local_Sex_Chrom_Dups': local_sex_count,
        'Intra_Chrom_Dups': intra_count,
        'Intra_Sex_Chrom_Dups': intra_sex_count,
        'Inter_Chrom_Dups': inter_chrom_count,
        'Inter_Sex_Chrom_Dups': inter_sex_count,
        'Conserved_Inter_Chrom_Genes': conserved_wgd_str,
        'Sex_Linked_Dups': sex_linked_count,
        'Conserved_Sex_Linked_Genes': conserved_sex_str,
        'Potential_Inter_Chrom_Dups' : pot_inter_chrom_count,
        'Unclassified_Complex' : unclassified_count,
    })

# ---------------------------------------------------------
# 4. Export Summary
# ---------------------------------------------------------
summary_df = pd.DataFrame(summary_records)
# Sort by highest number of samples first
summary_df = summary_df.sort_values(by='Total_Samples_With_Dups', ascending=False)
summary_df.to_csv(out_summary, sep='\t', index=False)

print(f"Successfully summarized all duplication types for {len(summary_df)} Orders.")