import pandas as pd
import argparse
import os

def get_best_score(file_path, min_cov, min_id):
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None, 0, 0
    
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]
    df = pd.read_csv(file_path, sep='\t', names=cols)
    df['coverage'] = df['length'] / df['qlen']
    
    filtered = df[df['coverage'] >= min_cov]
    
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
    parser.add_argument("--min_id", type=float)
    args = parser.parse_args()

    m_contig, m_score, m_cov = get_best_score(args.male, args.min_cov, args.min_id)
    f_contig, f_score, f_cov = get_best_score(args.female, args.min_cov, args.min_id)

    # Sex identification logic: Select higher bit-score (male vs female)
    if m_score > f_score:
        res = {'contig_id': m_contig, 'sex': 'V', 'score': m_score, 'cov': m_cov}
    elif f_score > m_score:
        res = {'contig_id': f_contig, 'sex': 'U', 'score': f_score, 'cov': f_cov}
    else:
        res = {'contig_id': 'None', 'sex': 'Unknown', 'score': 0, 'cov': 0}

    pd.DataFrame([res]).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()