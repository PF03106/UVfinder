import pandas as pd
import argparse
import os
import re

def parse_chrom_id(header_string):
    """
    Extract chromosome identifier (number or string) coming after "Chr".
    e.g.: '...Chr1...' -> '1', '...ChrU...' -> 'U'
    If no match found, returns the original string.
    """
    
    match = re.search(r'Chr(\d{1,2}|[a-zA-Z])$', str(header_string), re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return str(header_string)

def load_sex_map(file_path):
    """
    Load sex map TSV file into a dictionary { Chrom_ID : Sex_Tag }.
    Sex Tag Mapping:
      - U -> 'U'
      - V -> 'V'
      - Unknown -> 'N'
      - Everything else -> 'A' (Autosome)
    """
    mapping = {}
    default_tag = 'A' # Default to Autosome
    
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return mapping, 'N' # If no file, everything is Unknown (N)
        
    try:
        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            sex_chromosome = str(row['sex_chromosome']).strip()
            sex = str(row['sex']).strip()
            
            # If sex is unknown or no sex chromosome was found, set default tag to 'N'
            if sex_chromosome == 'Unknown' or sex == 'Unknown':
                default_tag = 'N'
                continue 
            
            # Valid sex chromosome found
            c_id = parse_chrom_id(sex_chromosome)
            if sex == 'U': 
                mapping[c_id] = 'U'
            elif sex == 'V': 
                mapping[c_id] = 'V'
            
    except Exception as e:
        print(f"Warning loading sex map: {e}")
        default_tag = 'N' # Fallback to N on error
        
    return mapping, default_tag

def is_overlapping(start1, end1, start2, end2, threshold=0.5):
    """
    Check if two genomic regions overlap by a certain ratio.
    threshold: 0.5 means if they overlap by more than 50% of the shorter length, it's a duplicate.
    """
    # Normalize coordinates (ensure min is first) because BLAST hits can be reverse strand
    s1, e1 = sorted([int(start1), int(end1)])
    s2, e2 = sorted([int(start2), int(end2)])
    
    # Calculate intersection
    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)
    
    # If no overlap exists
    if overlap_end <= overlap_start:
        return False
    
    overlap_len = overlap_end - overlap_start
    len1 = e1 - s1
    len2 = e2 - s2
    
    # Check if overlap ratio exceeds threshold for either sequence
    if len1 > 0 and (overlap_len / len1 > threshold): return True
    if len2 > 0 and (overlap_len / len2 > threshold): return True
    
    return False

def partition_blast(blast_file, sex_map_file, out_best, out_all):
    # Create output directories if they don't exist
    os.makedirs(os.path.dirname(out_best), exist_ok=True)
    os.makedirs(os.path.dirname(out_all), exist_ok=True)

    # Standard BLAST output format 6 columns
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    # Handle empty BLAST file
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        empty_df = pd.DataFrame(columns=cols+['rank', 'sex_tag'])
        empty_df.to_csv(out_best, sep='\t', index=False)
        empty_df.to_csv(out_all, sep='\t', index=False)
        return

    # 1. Load BLAST Results
    df = pd.read_csv(blast_file, sep='\t', names=cols, usecols=range(12))
    
    # 2. Pre-processing & Sorting
    # Extract clean Locus ID (remove pipes e.g., G4471|... -> G4471)
    df['locus_temp'] = df['qseqid'].apply(lambda x: str(x).split('|')[0])
    
    # Sort: Locus (A-Z) -> Evalue (Low to High) -> Bitscore (High to Low)
    # This ensures the "Best" hit is always the first row for each locus.
    df = df.sort_values(by=['locus_temp', 'evalue', 'bitscore'], ascending=[True, True, False])
    
    # 3. Load Sex Map
    sex_map = load_sex_map(sex_map_file)
    
    best_rows = []
    all_rows = [] # This will store only non-redundant hits
    
    # 4. Group by Locus and Filter
    # sort=False maintains the sorting order we just applied
    for locus, group in df.groupby('locus_temp', sort=False):
        
        # --- Process Best Hit ---
        # Since we sorted, the first row is the Best Hit
        best_hit = group.iloc[0].to_dict()
        best_chrom = parse_chrom_id(best_hit['sseqid'])
        best_hit['sex_tag'] = sex_map.get(best_chrom, 'A')
        best_hit['rank'] = 1
        
        best_rows.append(best_hit)
        
        # --- Process All Hits (with Overlap Filtering) ---
        saved_regions = [] # To track regions we have already accepted
        current_rank = 1
        
        for _, row in group.iterrows():
            chrom = row['sseqid']
            sstart, send = int(row['sstart']), int(row['send'])
            
            # Check for redundancy
            is_redundant = False
            for saved in saved_regions:
                saved_chrom, saved_start, saved_end = saved
                
                # If same chromosome AND overlaps significantly -> It is redundant
                if chrom == saved_chrom and is_overlapping(sstart, send, saved_start, saved_end):
                    is_redundant = True
                    break
            
            # If it's a new, distinct location (Real Duplication or Best Hit)
            if not is_redundant:
                saved_regions.append((chrom, sstart, send))
                
                row_dict = row.to_dict()
                chrom_id = parse_chrom_id(chrom)
                row_dict['sex_tag'] = sex_map.get(chrom_id, 'A')
                row_dict['rank'] = current_rank
                
                # Clean up temp column before saving
                if 'locus_temp' in row_dict: del row_dict['locus_temp']
                
                all_rows.append(row_dict)
                current_rank += 1
            
    # 5. Save Results
    out_cols = cols + ['rank', 'sex_tag']
    
    # Save Best Hits
    if best_rows:
        pd.DataFrame(best_rows)[out_cols].to_csv(out_best, sep='\t', index=False)
    else:
        pd.DataFrame(columns=out_cols).to_csv(out_best, sep='\t', index=False)
        
    # Save All Hits (Filtered)
    if all_rows:
        pd.DataFrame(all_rows)[out_cols].to_csv(out_all, sep='\t', index=False)
    else:
        pd.DataFrame(columns=out_cols).to_csv(out_all, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True, help="Input BLAST result file")
    parser.add_argument("--sex_map", required=True, help="Input Sex assignment TSV file")
    parser.add_argument("--out_best", required=True, help="Output path for Best Hits TSV")
    parser.add_argument("--out_all", required=True, help="Output path for All Hits TSV")
    args = parser.parse_args()
    
    partition_blast(args.blast, args.sex_map, args.out_best, args.out_all)