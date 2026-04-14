#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def get_best_score(file_path, min_bitscore_ratio, marker_total_length):
    """
    Parses BLAST output, filters noise, sums bitscores by contig, 
    and returns the best contig's score normalized by the actual marker length.
    """
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return "None", 0, marker_total_length
    
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    try:
        df = pd.read_csv(file_path, sep='\t', names=cols, dtype={'sseqid': str})
    except Exception:
        return "None", 0, marker_total_length

    if df.empty:
        return "None", 0, marker_total_length

    # Step 1: Pre-summation Noise Filter (Remove weak random hits)
    df = df[df['bitscore'] >= 50]
    
    if df.empty:
        return "None", 0, marker_total_length

    # Step 2: Sum bitscores by contig (Clustering multiple exons)
    grouped = df.groupby('sseqid')['bitscore'].sum().reset_index()
    
    # Step 3: Post-summation Ratio Filter
    absolute_max = grouped['bitscore'].max()
    cutoff = absolute_max * min_bitscore_ratio
    filtered = grouped[grouped['bitscore'] >= cutoff]
    
    if filtered.empty:
        return "None", 0, marker_total_length
    
    # Step 4: Select the contig with the highest summed score
    best = filtered.sort_values(by='bitscore', ascending=False).iloc[0]
    
    return best['sseqid'], best['bitscore'], marker_total_length


def main():
    parser = argparse.ArgumentParser(description="Sex ID based on clustered & length-normalized BLAST bitscores.")
    parser.add_argument("--sample", required=True, help="Target Sample ID")
    parser.add_argument("--samples_tsv", required=True, help="Metadata TSV")
    parser.add_argument("--male_blast", required=True, help="Male BLAST output (outfmt 6)")
    parser.add_argument("--female_blast", required=True, help="Female BLAST output (outfmt 6)")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--min_bitscore_ratio_UV", type=float, default=0.8, help="Threshold ratio (default: 0.8)")
    args = parser.parse_args()

    # Hardcoded actual marker lengths (in bp) for fair normalization
    FEMALE_MARKER_LEN = 586
    MALE_MARKER_LEN = 379

    # 1. Load metadata to check taxonomic Order
    try:
        meta_df = pd.read_csv(args.samples_tsv, sep='\t').set_index("sample_id")
        sample_order = str(meta_df.loc[args.sample, "order"]).lower().strip()
    except Exception:
        sample_order = "unknown"

    # 2. Get clustered scores and lengths
    m_contig, m_score, m_length = get_best_score(args.male_blast, args.min_bitscore_ratio_UV, MALE_MARKER_LEN)
    f_contig, f_score, f_length = get_best_score(args.female_blast, args.min_bitscore_ratio_UV, FEMALE_MARKER_LEN)

    # 3. Normalize scores by total marker length
    m_norm = (m_score / m_length) if m_length > 0 else 0
    f_norm = (f_score / f_length) if f_length > 0 else 0

    # 4. Initialize result dictionary
    res = {
        'sample_id': args.sample,
        'sex': 'Unknown',
        'sex_chromosome': 'Unknown',
        'sum_bit_score': 0,
        'basis': 'No_Signal'
    }

    # --- Classification Logic ---

    # Case A: Female signal is stronger (Normalized)
    if f_norm > m_norm and f_score > 0:
        res.update({
            'sex': 'U',
            'sex_chromosome': f_contig,
            'sum_bit_score': f_score,
            'basis': 'Confirmed_Female_Hit'
        })
    
    # Case B: Male signal is stronger (Normalized)
    elif m_norm > f_norm and m_score > 0:
        res.update({
            'sex': 'V',
            'sex_chromosome': m_contig,
            'sum_bit_score': m_score,
            'basis': 'Confirmed_Male_Hit'
        })
    
    # Case C: No signal from either marker
    elif m_score == 0 and f_score == 0:
        if "sphagnales" in sample_order:
            res.update({
                'sex': 'V',
                'sex_chromosome': 'Unknown',
                'sum_bit_score': 0,
                'basis': 'Inferred_by_Order_Sphagnales'
            })
        else:
            res.update({
                'sex': 'Unknown',
                'sex_chromosome': 'Unknown',
                'sum_bit_score': 0,
                'basis': 'No_Signal_Non_Sphagnales'
            })
    
    # Case D: Equal normalized scores (Ambiguous/Monoecious)
    else:
        res.update({
            'sex': 'Unknown',
            'sex_chromosome': 'Ambiguous',
            'sum_bit_score': max(m_score, f_score),
            'basis': 'Ambiguous_Equal_Scores'
        })

    # 5. Save the output to TSV
    pd.DataFrame([res]).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()