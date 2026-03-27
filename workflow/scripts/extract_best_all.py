import pandas as pd
import argparse
import os
import re

def parse_chrom_id(header_string):
    """
    Extract the chromosome identifier (digit 1-2 or single letter) after 'Chr'.
    e.g., 'scaffold_Chr1' -> '1', 'ChrU' -> 'U'
    Uses a capturing group (\d{1,2}|[a-zA-Z]) to avoid IndexError.
    """
    # Pattern looks for 'Chr' followed by 1-2 digits OR a single letter at the end of the string
    match = re.search(r'Chr(\d{1,2}|[a-zA-Z])$', str(header_string), re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return str(header_string)

def load_sex_map(file_path):
    """
    Load the sex assignment TSV and create a mapping dictionary.
    Returns: (dict of {Chrom_ID: Tag}, default_tag_string)
    """
    mapping = {}
    default_tag = 'A' # 'A' for Autosome
    
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return mapping, 'N' # 'N' for No information/Unknown
        
    try:
        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            sex_chromosome = str(row['sex_chromosome']).strip()
            sex = str(row['sex']).strip()
            
            # If the sample is labeled 'Unknown', set the default tag to 'N' for all hits
            if sex_chromosome == 'Unknown' or sex == 'Unknown':
                default_tag = 'N'
                continue 
            
            # Map specific chromosome IDs to their sex tags (U or V)
            c_id = parse_chrom_id(sex_chromosome)
            if sex == 'U': 
                mapping[c_id] = 'U'
            elif sex == 'V': 
                mapping[c_id] = 'V'
            
    except Exception as e:
        print(f"Warning loading sex map: {e}")
        default_tag = 'N'
        
    return mapping, default_tag

def is_overlapping(start1, end1, start2, end2, threshold=0.5):
    """
    Determine if two genomic intervals overlap significantly.
    Used to filter redundant BLAST hits on the same locus.
    """
    s1, e1 = sorted([int(start1), int(end1)])
    s2, e2 = sorted([int(start2), int(end2)])
    
    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)
    
    if overlap_end <= overlap_start:
        return False
    
    overlap_len = overlap_end - overlap_start
    len1 = e1 - s1
    len2 = e2 - s2
    
    # Check if overlap exceeds the threshold for either hit
    if len1 > 0 and (overlap_len / len1 > threshold): return True
    if len2 > 0 and (overlap_len / len2 > threshold): return True
    
    return False

def partition_blast(blast_file, sex_map_file, out_best, out_all, overlap_threshold):
    """
    Main logic to tag BLAST hits and partition them into Best and All hits.
    """
    os.makedirs(os.path.dirname(out_best), exist_ok=True)
    os.makedirs(os.path.dirname(out_all), exist_ok=True)

    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        empty_df = pd.DataFrame(columns=cols+['rank', 'sex_tag'])
        empty_df.to_csv(out_best, sep='\t', index=False)
        empty_df.to_csv(out_all, sep='\t', index=False)
        return

    # 1. Load BLAST Results
    df = pd.read_csv(blast_file, sep='\t', names=cols, usecols=range(12))
    
    # 2. Extract clean Locus ID and sort by quality
    df['locus_temp'] = df['qseqid'].apply(lambda x: str(x).split('|')[0])
    df = df.sort_values(by=['locus_temp', 'evalue', 'bitscore'], ascending=[True, True, False])
    
    # 3. Load Sex Map (Unpack the tuple)
    sex_map, default_tag = load_sex_map(sex_map_file) 
    
    best_rows = []
    all_rows = [] 
    
    # 4. Process each locus
    for locus, group in df.groupby('locus_temp', sort=False):
        
        # --- Handle Best Hit ---
        best_hit = group.iloc[0].to_dict()
        best_chrom = parse_chrom_id(best_hit['sseqid'])
        best_hit['sex_tag'] = sex_map.get(best_chrom, default_tag)
        best_hit['rank'] = 1
        best_rows.append(best_hit)
        
        # --- Handle All Hits (Filtering Redundant Overlaps) ---
        saved_regions = [] 
        current_rank = 1
        
        for _, row in group.iterrows():
            chrom = row['sseqid']
            sstart, send = int(row['sstart']), int(row['send'])
            
            is_redundant = False
            for saved in saved_regions:
                saved_chrom, saved_start, saved_end = saved
                if chrom == saved_chrom and is_overlapping(sstart, send, saved_start, saved_end, threshold=overlap_threshold):
                    is_redundant = True
                    break
            
            if not is_redundant:
                saved_regions.append((chrom, sstart, send))
                row_dict = row.to_dict()
                chrom_id = parse_chrom_id(chrom)
                
                # Tag hit: either U/V from map or default_tag (A/N)
                row_dict['sex_tag'] = sex_map.get(chrom_id, default_tag)
                row_dict['rank'] = current_rank
                
                if 'locus_temp' in row_dict: del row_dict['locus_temp']
                all_rows.append(row_dict)
                current_rank += 1
            
    # 5. Save Results
    out_cols = cols + ['rank', 'sex_tag']
    pd.DataFrame(best_rows)[out_cols].to_csv(out_best, sep='\t', index=False) if best_rows else pd.DataFrame(columns=out_cols).to_csv(out_best, sep='\t', index=False)
    pd.DataFrame(all_rows)[out_cols].to_csv(out_all, sep='\t', index=False) if all_rows else pd.DataFrame(columns=out_cols).to_csv(out_all, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tag BLAST hits with sex information and filter overlaps.")
    parser.add_argument("--blast", required=True, help="Input BLAST result file (tsv)")
    parser.add_argument("--sex_map", required=True, help="Input sex assignment file (tsv)")
    parser.add_argument("--out_best", required=True, help="Output path for best hits per locus")
    parser.add_argument("--out_all", required=True, help="Output path for all non-redundant hits")
    parser.add_argument("--overlap_threshold", type=float, default=0.5, help="Threshold for overlapping BLAST hits (0.0 to 1.0)")
    args = parser.parse_args()
    
    partition_blast(args.blast, args.sex_map, args.out_best, args.out_all, overlap_threshold=args.overlap_threshold)
