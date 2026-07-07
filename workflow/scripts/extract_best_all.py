import pandas as pd
import argparse
import os
import re

def parse_chrom_id(header_string):
    """
    Extract the chromosome identifier (digit 1-2 or single letter) after 'Chr'.
    e.g., 'scaffold_Chr1' -> '1', 'ChrU' -> 'U'
    """
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
            
            # If the sample is labeled 'Unknown', set the default tag to 'N'
            if sex_chromosome == 'Unknown' or sex == 'Unknown':
                default_tag = 'N'
                continue 
            
            # Map specific chromosome IDs to their sex tags
            c_id = parse_chrom_id(sex_chromosome)
            if sex == 'U': mapping[c_id] = 'U'
            elif sex == 'V': mapping[c_id] = 'V'
                
    except Exception as e:
        print(f"Warning loading sex map: {e}")
        default_tag = 'N'
        
    return mapping, default_tag

def is_overlapping(start1, end1, start2, end2, threshold=0.5):
    """
    Determine if two genomic intervals overlap significantly based on a threshold.
    """
    s1, e1 = sorted([int(start1), int(end1)])
    s2, e2 = sorted([int(start2), int(end2)])
    
    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)
    
    if overlap_end <= overlap_start: return False
    
    overlap_len = overlap_end - overlap_start
    len1 = e1 - s1
    len2 = e2 - s2
    
    # Check if overlap exceeds the threshold for either hit
    if len1 > 0 and (overlap_len / len1 > threshold): return True
    if len2 > 0 and (overlap_len / len2 > threshold): return True
    return False

def partition_blast(blast_file, sex_map_file, out_best, out_all, overlap_threshold, min_pident, min_bitscore_ratio):
    os.makedirs(os.path.dirname(out_best), exist_ok=True)
    os.makedirs(os.path.dirname(out_all), exist_ok=True)

    # Maintain exactly 14 columns for the output
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    out_cols = cols + ['rank', 'sex_tag']
    
    def save_empty():
        pd.DataFrame(columns=out_cols).to_csv(out_best, sep='\t', index=False)
        pd.DataFrame(columns=out_cols).to_csv(out_all, sep='\t', index=False)

    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        save_empty()
        return

    df = pd.read_csv(blast_file, sep='\t', names=cols, usecols=range(12))
    if df.empty:
        save_empty()
        return

    # Force numeric conversion for coordinates to prevent errors
    df['sstart'] = pd.to_numeric(df['sstart'], errors='coerce')
    df['send'] = pd.to_numeric(df['send'], errors='coerce')
    df = df.dropna(subset=['sstart', 'send'])

    # Extract Locus ID from qseqid
    df['locus_temp'] = df['qseqid'].apply(lambda x: str(x).split('|')[0])

    # 1. Filter by percent identity (pident)
    df = df[df['pident'] >= min_pident]

    if not df.empty:
        # 2. Filter by bitscore ratio (relative to the maximum bitscore of that locus)
        max_bitscores = df.groupby('locus_temp')['bitscore'].transform('max')
        df = df[df['bitscore'] >= (max_bitscores * min_bitscore_ratio)]

    if df.empty:
        save_empty()
        return

    # Sort by highest quality: evalue ascending (lower is better), bitscore descending (higher is better)
    df = df.sort_values(by=['locus_temp', 'evalue', 'bitscore'], ascending=[True, True, False])
    
    sex_map, default_tag = load_sex_map(sex_map_file) 
    
    best_rows = []
    all_rows = [] 
    
    for locus, group in df.groupby('locus_temp', sort=False):
        
        # --- Handle Best Hit ---
        best_hit = group.iloc[0].to_dict()
        best_chrom = parse_chrom_id(best_hit['sseqid'])
        best_hit['sex_tag'] = sex_map.get(best_chrom, default_tag)
        best_hit['rank'] = 1
        if 'locus_temp' in best_hit: del best_hit['locus_temp']
        best_rows.append(best_hit)
        
        # --- Handle All Hits (Collect valid non-overlapping hits) ---
        saved_regions = [] 
        valid_hits = []
        current_rank = 1
        
        for _, row in group.iterrows():
            chrom = row['sseqid']
            sstart, send = int(row['sstart']), int(row['send'])
            
            # Check for physical overlap with already saved regions
            is_redundant = False
            for saved in saved_regions:
                saved_chrom, saved_start, saved_end = saved
                if chrom == saved_chrom and is_overlapping(sstart, send, saved_start, saved_end, threshold=overlap_threshold):
                    is_redundant = True
                    break
            
            # Save if it is a non-overlapping valid hit
            if not is_redundant:
                saved_regions.append((chrom, sstart, send))
                row_dict = row.to_dict()
                chrom_id = parse_chrom_id(chrom)
                
                row_dict['sex_tag'] = sex_map.get(chrom_id, default_tag)
                row_dict['rank'] = current_rank
                if 'locus_temp' in row_dict: del row_dict['locus_temp']
                
                valid_hits.append(row_dict)
                current_rank += 1
                
        # --- Duplication Determination and Storage Logic ---
        if not valid_hits:
            continue
            
        if len(valid_hits) == 1:
            # If there's only one valid hit, append it to all_rows
            all_rows.append(valid_hits[0])
        else:
            # Calculate frequency of each probe name (qseqid) among valid hits
            qseqid_counts = {}
            for h in valid_hits:
                q = h['qseqid']
                qseqid_counts[q] = qseqid_counts.get(q, 0) + 1
                
            # Identify true duplication probes appearing 2 or more times
            dup_qseqids = set([q for q, count in qseqid_counts.items() if count > 1])
            
            if dup_qseqids:
                # Keep all hits for these duplicated probes, and ALWAYS ensure the best hit (rank 1) is included
                for h in valid_hits:
                    if h['qseqid'] in dup_qseqids or h['rank'] == 1:
                        all_rows.append(h)
            else:
                # If all probes hit exactly once (false duplication or different plants), keep only the best hit
                all_rows.append(valid_hits[0])

    # Save Results to TSV format
    if best_rows:
        pd.DataFrame(best_rows)[out_cols].to_csv(out_best, sep='\t', index=False)
    else:
        pd.DataFrame(columns=out_cols).to_csv(out_best, sep='\t', index=False)
        
    if all_rows:
        pd.DataFrame(all_rows)[out_cols].to_csv(out_all, sep='\t', index=False)
    else:
        pd.DataFrame(columns=out_cols).to_csv(out_all, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tag BLAST hits, filter, get best hits and specific duplicated hits.")
    parser.add_argument("--blast", required=True, help="Input BLAST result file (tsv)")
    parser.add_argument("--sex_map", required=True, help="Input sex assignment file (tsv)")
    parser.add_argument("--out_best", required=True, help="Output path for best hits per locus")
    parser.add_argument("--out_all", required=True, help="Output path for all non-redundant/duplicated hits")
    parser.add_argument("--overlap_threshold", type=float, default=0.5, help="Threshold for overlapping BLAST hits (0.0 to 1.0)")
    parser.add_argument("--min_pident", type=float, default=70.0, help="Minimum percent identity to keep a hit")
    parser.add_argument("--min_bitscore_ratio", type=float, default=0.8, help="Minimum bitscore ratio")
    parser.add_argument("--max_dist", type=int, default=100000, help="Ignored in this script, kept for pipeline compatibility") 
    
    args = parser.parse_args()
    
    partition_blast(
        blast_file=args.blast, 
        sex_map_file=args.sex_map, 
        out_best=args.out_best, 
        out_all=args.out_all, 
        overlap_threshold=args.overlap_threshold,
        min_pident=args.min_pident,
        min_bitscore_ratio=args.min_bitscore_ratio
    )