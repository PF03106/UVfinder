#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def get_best_score(file_path, min_bitscore_ratio):
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None, 0
    
    # Load columns without qlen since coverage is not calculated
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    try:
        df = pd.read_csv(file_path, sep='\t', names=cols)
    except pd.errors.EmptyDataError:
        return None, 0

    if df.empty:
        return None, 0

    # Calculate the absolute maximum bitscore in the entire file
    absolute_max = df['bitscore'].max()

    # 1. Filter hits that are >= N% (e.g., 80%) of the absolute maximum bitscore
    filtered = df[df['bitscore'] >= (absolute_max * min_bitscore_ratio)]
    
    # 2. Keep only hits that match the formal Chromosome (Chr) naming convention
    filtered = filtered[filtered['sseqid'].str.contains(r'Chr(\d{1,2}|[a-zA-Z])$', regex=True, na=False)]
    
    if filtered.empty:
        return None, 0
    
    # 3. Return the top hit with the highest bitscore among the filtered Chr hits
    best = filtered.sort_values(by='bitscore', ascending=False).iloc[0]
    return best['sseqid'], best['bitscore']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--male", required=True)
    parser.add_argument("--female", required=True)
    parser.add_argument("--output", required=True)
    # Bitscore ratio threshold (default 0.8)
    parser.add_argument("--min_bitscore_ratio", type=float, default=0.8)
    args = parser.parse_args()

    m_contig, m_score = get_best_score(args.male, args.min_bitscore_ratio)
    f_contig, f_score = get_best_score(args.female, args.min_bitscore_ratio)

    # Sex identification logic (Determine U/V)
    if m_score > f_score:
        res = {'sex_chromosome': m_contig, 'sex': 'V', 'bit_score': m_score}
    elif f_score > m_score:
        res = {'sex_chromosome': f_contig, 'sex': 'U', 'bit_score': f_score}
    else:
        res = {'sex_chromosome': 'Unknown', 'sex': 'Unknown', 'bit_score': 0}

    pd.DataFrame([res]).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()