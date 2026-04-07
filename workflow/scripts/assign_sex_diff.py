#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def get_best_score(file_path, min_bitscore_ratio_UV):
    """
    Parses BLAST output to find the best bitscore among hits that pass 
    the bitscore ratio threshold and match chromosome naming conventions.
    """
    # Check if file exists and is not empty
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return "None", 0
    
    # BLAST outfmt 6 columns
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    try:
        df = pd.read_csv(file_path, sep='\t', names=cols, dtype={'sseqid': str})
    except Exception:
        return "None", 0

    if df.empty:
        return "None", 0

    # Calculate filtering threshold based on the top hit in the file
    absolute_max = df['bitscore'].max()
    cutoff = absolute_max * min_bitscore_ratio_UV
    
    # Step 1: Filter hits based on the bitscore ratio
    filtered = df[df['bitscore'] >= cutoff]
    
    # Step 2: Filter for formal Chromosome (Chr) naming convention
    # Adjust regex if your assembly uses different naming (e.g., 'scaffold')
    filtered = filtered[filtered['sseqid'].str.contains(r'Chr(\d{1,2}|[a-zA-Z])$', regex=True, na=False)]
    
    if filtered.empty:
        return "None", 0
    
    # Step 3: Select the hit with the highest bitscore
    best = filtered.sort_values(by='bitscore', ascending=False).iloc[0]
    return best['sseqid'], best['bitscore']

def main():
    parser = argparse.ArgumentParser(description="Sex identification based on normalized BLAST bitscores.")
    parser.add_argument("--sample", required=True, help="Target Sample ID")
    parser.add_argument("--samples_tsv", required=True, help="Metadata TSV containing 'order' information")
    parser.add_argument("--male_blast", required=True, help="BLAST result for Male probes")
    parser.add_argument("--female_blast", required=True, help="BLAST result for Female probes")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    parser.add_argument("--min_bitscore_ratio_UV", type=float, default=0.8, help="Filtering threshold (default: 0.8)")
    args = parser.parse_args()

    # 1. Load metadata to check taxonomic Order
    try:
        meta_df = pd.read_csv(args.samples_tsv, sep='\t').set_index("sample_id")
        sample_order = str(meta_df.loc[args.sample, "order"]).lower().strip()
    except Exception as e:
        print(f"Error loading metadata for {args.sample}: {e}")
        sample_order = "unknown"

    # 2. Retrieve best bitscores from BLAST results
    m_contig, m_score = get_best_score(args.male_blast, args.min_bitscore_ratio_UV)
    f_contig, f_score = get_best_score(args.female_blast, args.min_bitscore_ratio_UV)

    # 3. Normalize scores by marker length
    # Male marker: 384 bp | Female marker: 595 bp
    m_norm = m_score / 378
    f_norm = f_score / 586

    # 4. Initialize result dictionary
    res = {
        'sample_id': args.sample,
        'sex': 'Unknown',
        'sex_chromosome': 'Unknown',
        'bit_score': 0,
        'basis': 'No_Signal'
    }

    # --- Classification Logic (Rule-based) ---

    # Case A: Female signal is stronger (Normalized)
    if f_norm > m_norm and f_score > 0:
        res.update({
            'sex': 'U',
            'sex_chromosome': f_contig,
            'bit_score': f_score,
            'basis': 'Confirmed_Female_Hit'
        })
    
    # Case B: Male signal is stronger (Normalized)
    elif m_norm > f_norm and m_score > 0:
        res.update({
            'sex': 'V',
            'sex_chromosome': m_contig,
            'bit_score': m_score,
            'basis': 'Confirmed_Male_Hit'
        })
    
    # Case C: No signal from either marker (Handling Sphagnales exception)
    elif m_score == 0 and f_score == 0:
        if "sphagnales" in sample_order:
            # Special inference for Sphagnales: absence of female signal usually implies Male (V)
            res.update({
                'sex': 'V',
                'sex_chromosome': 'Unknown',
                'bit_score': 0,
                'basis': 'Inferred_by_Order_Sphagnales'
            })
        else:
            # Generic case for other taxa with no signal
            res.update({
                'sex': 'Unknown',
                'sex_chromosome': 'Unknown',
                'bit_score': 0,
                'basis': 'No_Signal_Non_Sphagnales'
            })
    
    # Case D: Equal normalized scores (very rare)
    else:
        res.update({
            'sex': 'Unknown',
            'sex_chromosome': 'Ambiguous',
            'bit_score': max(m_score, f_score),
            'basis': 'Ambiguous_Equal_Scores'
        })

    # 5. Save the identification result to a TSV file
    pd.DataFrame([res]).to_csv(args.output, sep='\t', index=False)
    print(f"Done: {args.sample} identified as {res['sex']} via {res['basis']}")

if __name__ == "__main__":
    main()