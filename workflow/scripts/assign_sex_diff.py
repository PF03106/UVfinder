import pandas as pd
import argparse
import os

def get_best_score(file_path, min_cov, max_evalue):
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None, 0, 0
    
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]
    df = pd.read_csv(file_path, sep='\t', names=cols)
    df['coverage'] = df['length'] / df['qlen']
    
    filtered = df[(df['coverage'] >= min_cov) & (df['evalue'] <= max_evalue)]
    
    if filtered.empty:
        return None, 0, 0
    
    # Return Top Hit which has the highest Bit-score
    best = filtered.sort_values(by='bitscore', ascending=False).iloc[0]
    return best['sseqid'], best['bitscore'], best['coverage']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--male", required=True)
    parser.add_argument("--female", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min_cov", type=float)
    parser.add_argument("--max_evalue", type=float)
    args = parser.parse_args()

    m_contig, m_score, m_cov = get_best_score(args.male, args.min_cov, args.max_evalue)
    f_contig, f_score, f_cov = get_best_score(args.female, args.min_cov, args.max_evalue)

    # Sex identification logic: Select higher bit-score (male vs female) (m_score and f_score are bit-scores of male and female best hits respectively)
    if m_score > f_score:
        res = {'sex_chromosome': m_contig, 'sex': 'V', 'bit_score': m_score, 'coverage': m_cov}
    elif f_score > m_score:
        res = {'sex_chromosome': f_contig, 'sex': 'U', 'bit_score': f_score, 'coverage': f_cov}
    else:
        res = {'sex_chromosome': 'Unknown', 'sex': 'Unknown', 'bit_score': 0, 'coverage': 0}

    pd.DataFrame([res]).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()