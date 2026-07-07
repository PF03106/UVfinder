#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def get_fasta_lengths(fasta_path):
    """
    Reads a FASTA file and returns a dictionary mapping sequence IDs to their lengths.
    """
    lengths = {}
    if not os.path.exists(fasta_path):
        return lengths
        
    with open(fasta_path, 'r') as f:
        seq_id = None
        length = 0
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id is not None:
                    lengths[seq_id] = length
                # Get the ID (first word after '>')
                seq_id = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
        if seq_id is not None:
            lengths[seq_id] = length
            
    return lengths

def get_best_score(file_path, min_bitscore_ratio, marker_lengths):
    """
    Parses BLAST output, groups by BOTH contig and marker ID to find the best 
    single marker hit, and returns the contig, marker ID, score, and specific marker length.
    """
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return "None", "None", 0, 1 # default length 1 to prevent division by zero
    
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    try:
        df = pd.read_csv(file_path, sep='\t', names=cols, dtype={'sseqid': str, 'qseqid': str})
    except Exception:
        return "None", "None", 0, 1

    if df.empty:
        return "None", "None", 0, 1

    # Step 1: Pre-summation Noise Filter (Remove weak random hits)
    df = df[df['bitscore'] >= 50]
    
    if df.empty:
        return "None", "None", 0, 1

    # Step 2: Sum bitscores by both contig AND specific marker (qseqid)
    # This ensures if multiple markers were used, we evaluate them independently
    grouped = df.groupby(['sseqid', 'qseqid'])['bitscore'].sum().reset_index()
    
    # Step 3: Post-summation Ratio Filter
    absolute_max = grouped['bitscore'].max()
    cutoff = absolute_max * min_bitscore_ratio
    filtered = grouped[grouped['bitscore'] >= cutoff]
    
    if filtered.empty:
        return "None", "None", 0, 1
    
    # Step 4: Select the (contig, marker) pair with the highest summed score
    best = filtered.sort_values(by='bitscore', ascending=False).iloc[0]
    
    best_sseqid = best['sseqid']
    best_qseqid = best['qseqid']
    best_score = best['bitscore']
    
    # Get the exact length for this specific winning marker
    best_marker_length = marker_lengths.get(best_qseqid, 1) 
    
    return best_sseqid, best_qseqid, best_score, best_marker_length


def main():
    parser = argparse.ArgumentParser(description="Sex ID based on clustered & length-normalized BLAST bitscores.")
    parser.add_argument("--sample", required=True, help="Target Sample ID")
    parser.add_argument("--samples_tsv", required=True, help="Metadata TSV")
    parser.add_argument("--male_blast", required=True, help="Male BLAST output (outfmt 6)")
    parser.add_argument("--female_blast", required=True, help="Female BLAST output (outfmt 6)")
    parser.add_argument("--male_marker", required=True, help="Male marker FASTA file")
    parser.add_argument("--female_marker", required=True, help="Female marker FASTA file")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--min_bitscore_ratio_UV", type=float, default=0.8, help="Threshold ratio (default: 0.8)")
    args = parser.parse_args()

    # 1. Parse marker lengths dynamically from FASTA files
    f_lengths_dict = get_fasta_lengths(args.female_marker)
    m_lengths_dict = get_fasta_lengths(args.male_marker)

    # 2. Load metadata to check taxonomic Order
    try:
        meta_df = pd.read_csv(args.samples_tsv, sep='\t').set_index("sample_id")
        sample_order = str(meta_df.loc[args.sample, "order"]).lower().strip()
    except Exception:
        sample_order = "unknown"

    # 3. Get best clustered scores, chosen marker ID, and specific lengths
    m_contig, m_marker, m_score, m_length = get_best_score(args.male_blast, args.min_bitscore_ratio_UV, m_lengths_dict)
    f_contig, f_marker, f_score, f_length = get_best_score(args.female_blast, args.min_bitscore_ratio_UV, f_lengths_dict)

    # 4. Normalize scores by the chosen marker's actual length
    m_norm = (m_score / m_length) if m_score > 0 else 0
    f_norm = (f_score / f_length) if f_score > 0 else 0

    # 5. Initialize result dictionary (added best_marker for tracking)
    res = {
        'sample_id': args.sample,
        'sex': 'Unknown',
        'sex_chromosome': 'Unknown',
        'sum_bit_score': 0,
        'basis': 'No_Signal',
        'best_female_marker': f_marker if f_marker != "None" else "NA",
        'best_male_marker': m_marker if m_marker != "None" else "NA"
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

    # 6. Save the output to TSV
    pd.DataFrame([res]).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()