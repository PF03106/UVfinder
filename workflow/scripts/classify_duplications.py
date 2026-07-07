import pandas as pd
import sys
import os
import argparse

# ---------------------------------------------------------
# 0. Argument Parsing 
# ---------------------------------------------------------
parser = argparse.ArgumentParser(description="Classify multi-hit probes into duplication types with deduplication.")
parser.add_argument("--input", required=True, help="Input BLAST hits TSV file (filtered 'All' hits)")
parser.add_argument("--output", required=True, help="Output classified TSV file")
parser.add_argument("--max_dist", type=int, default=100000, help="Maximum distance for local duplication (bp)")
args = parser.parse_args()

input_file = args.input
output_file = args.output
max_dist = args.max_dist

# 1. Define BLAST output columns (outfmt 6 + rank + sex_tag)
cols = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'rank', 'sex_tag'
]

# Define output columns (keeping qseqid to identify the specific probe used)
out_cols = ['qseqid', 'duplication_type', 'hit_count', 'mapped_chromosomes', 'sex_tags']

# Safety Check 1: Check if input file is completely empty
if os.stat(input_file).st_size == 0:
    empty_df = pd.DataFrame(columns=out_cols)
    empty_df.to_csv(output_file, sep='\t', index=False)
    sys.exit(0)

# 2. Load data
df = pd.read_csv(input_file, sep='\t', names=cols)

# Force numeric conversion for coordinates to prevent string subtraction errors
df['sstart'] = pd.to_numeric(df['sstart'], errors='coerce')
df['send'] = pd.to_numeric(df['send'], errors='coerce')
df = df.dropna(subset=['sstart', 'send'])

# Safety Check 2: Check if dataframe is empty after loading/filtering
if df.empty:
    empty_df = pd.DataFrame(columns=out_cols)
    empty_df.to_csv(output_file, sep='\t', index=False)
    sys.exit(0)

# ---------------------------------------------------------
# 3. Filter for queries with multiple hits
# ---------------------------------------------------------
hit_counts = df.groupby('qseqid').size()
multi_queries = hit_counts[hit_counts > 1].index
df_multi = df[df['qseqid'].isin(multi_queries)]

results = []

if not df_multi.empty:
    # 4. Classification Logic
    for qseqid, group in df_multi.groupby('qseqid'):
        # Extract unique chromosomes and SORT them to ensure identical sets are matched
        # e.g., 'Chr1,Chr2' will be treated the same as 'Chr2,Chr1'
        unique_sseqids = sorted(list(group['sseqid'].unique()))
        unique_tags = set(group['sex_tag'].unique())
        
        # ---------------------------------------------------------
        # [Case 1] Hits mapped to the EXACT SAME chromosome/scaffold
        # ---------------------------------------------------------
        if len(unique_sseqids) == 1:
            min_pos = group[['sstart', 'send']].min().min()
            max_pos = group[['sstart', 'send']].max().max()
            distance = max_pos - min_pos
            
            has_sex = ('U' in unique_tags) or ('V' in unique_tags)
            
            if has_sex:
                dup_type = 'Local_Sex_Chrom_Duplication' if distance <= max_dist else 'Intra_Sex_Chrom_Duplication'
            else:
                dup_type = 'Local_Duplication' if distance <= max_dist else 'Intra_chromosomal_Duplication'
                
        # ---------------------------------------------------------
        # [Case 2] Hits mapped to MULTIPLE DIFFERENT chromosomes
        # ---------------------------------------------------------
        else:
            has_sex = ('U' in unique_tags) or ('V' in unique_tags)
            has_auto = 'A' in unique_tags
            has_unknown = 'N' in unique_tags
            
            if has_sex and has_auto:
                dup_type = 'Sex_Linked_Duplication'
            elif has_auto and not has_sex and not has_unknown:
                dup_type = 'Inter_Chromosomal_Duplication'
            elif has_sex and not has_auto and not has_unknown:
                dup_type = 'Inter_Sex_Chrom_Duplication'
            elif has_unknown:
                dup_type = 'Potential_Inter_Chrom_Dups'
            else:
                dup_type = 'Unclassified_Complex'

        # Append to results
        results.append({
            'qseqid': qseqid,
            'duplication_type': dup_type,
            'hit_count': len(group),
            'mapped_chromosomes': ",".join(unique_sseqids),
            'sex_tags': ",".join(unique_tags)
        })

# ---------------------------------------------------------
# 5. Summarizing Logic 
# Drop duplicates if 'locus', 'duplication_type', and 'mapped_chromosomes' are all identical.
# This case leave only the best(first) hits
# ---------------------------------------------------------
out_df = pd.DataFrame(results, columns=out_cols)

if not out_df.empty:
    # Extract the base locus ID (e.g., 'G4992' from 'G4992|hypnales_...')
    out_df['locus'] = out_df['qseqid'].apply(lambda x: str(x).split('|')[0])
    
    # Drop duplicates if 'locus', 'duplication_type', and 'mapped_chromosomes' are all identical.
    # keep='first' ensures we retain the representative probe (usually the highest scoring one).
    out_df = out_df.drop_duplicates(subset=['locus', 'duplication_type', 'mapped_chromosomes'], keep='first')
    
    # Remove the temporary 'locus' column to restore the original 5-column format
    out_df = out_df.drop(columns=['locus'])

# Export results
out_df.to_csv(output_file, sep='\t', index=False)