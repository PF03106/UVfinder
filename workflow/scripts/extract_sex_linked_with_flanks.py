#!/usr/bin/env python3
# Extract sequences from genome based on BLAST hits in TSV, with flanking regions and hit limits

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os

def extract_sequences(tsv_file, genome_file, loci_list_file, output_dir, sample_id, flank_bp, max_hits):
    print(f"--- Extracting sequences for {sample_id} ---")
    print(f"    (Parameters: Flank={flank_bp}bp, Max Hits={max_hits})")

    # Create output directory first to prevent Snakemake MissingOutputException
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load Target Loci List (Can be Sex-linked or Not-sex-linked depending on input)
    if not os.path.exists(loci_list_file) or os.path.getsize(loci_list_file) == 0:
        print("⚠️ Warning: Loci list file is missing or empty. Exiting safely.")
        return

    with open(loci_list_file) as f:
        target_loci = set(line.strip() for line in f if line.strip())

    # 2. Check TSV file
    if not os.path.exists(tsv_file) or os.path.getsize(tsv_file) == 0:
        print("⚠️ Warning: BLAST TSV file is missing or empty. Exiting safely.")
        return

    # 3. Load TSV first
    df = pd.read_csv(tsv_file, sep='\t')
    if df.empty:
        print("⚠️ Warning: BLAST TSV is empty (only headers). Exiting safely.")
        return
        
    df['clean_locus'] = df['qseqid'].apply(lambda x: str(x).split('|')[0])
    
    # Filter loci based on target_loci (Keep only loci present in the provided list)
    df = df[df['clean_locus'].isin(target_loci)]
    print(f"Target loci to extract: {len(set(df['clean_locus'].unique()))}")
    
    if df.empty:
        print("⚠️ Warning: No matching loci found after filtering. Exiting safely.")
        return

    extracted_count = 0

    # 4. Create indexed genome reader
    try:
        genome_index = SeqIO.index(genome_file, "fasta")
        
        # 5. Extract Sequence per Gene
        for locus, group in df.groupby('clean_locus'):
            records = []
            
            sorted_group = group.sort_values('rank')
            hits_processed = 0
            
            for _, row in sorted_group.iterrows():
                if hits_processed >= max_hits:
                    break
                    
                sseqid = str(row['sseqid'])
                if sseqid not in genome_index: 
                    continue
                
                seq_record = genome_index[sseqid]
                full_seq = seq_record.seq
                chrom_len = len(full_seq)
                start, end = int(row['sstart']), int(row['send'])
                rank = row['rank']
                sex_tag = row['sex_tag']
                
                # Apply flanking region
                if start < end: # Forward
                    ext_start = max(0, start - 1 - flank_bp)
                    ext_end = min(chrom_len, end + flank_bp)
                    seq_frag = full_seq[ext_start:ext_end]
                    strand_mark = "+"
                else: # Reverse
                    ext_start = max(0, end - 1 - flank_bp)
                    ext_end = min(chrom_len, start + flank_bp)
                    seq_frag = full_seq[ext_start:ext_end].reverse_complement()
                    strand_mark = "-"
                
                new_id = f"{sample_id}_{locus}_R{rank}_{sex_tag}"
                desc = f"orig:{sseqid}:{start}-{end}({strand_mark}) flank:{flank_bp}"
                
                records.append(SeqRecord(seq_frag, id=new_id, description=desc))
                hits_processed += 1
                
            if records:
                out_path = os.path.join(output_dir, f"{str(locus)}.fasta")
                with open(out_path, "w") as f:
                    SeqIO.write(records, f, "fasta")
                    extracted_count += 1
                    
    except Exception as e:
        print(f"❌ Error during genome indexing or extraction: {e}")
        return
    finally:
        if 'genome_index' in locals():
            genome_index.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--loci_list", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--sample", required=True)
    
    # Parameters to be received from Config
    parser.add_argument("--flank", type=int, default=20, help="Flanking bp")
    parser.add_argument("--max_hits", type=int, default=10, help="Max hits per gene")
    
    args = parser.parse_args()
    
    extract_sequences(args.tsv, args.genome, args.loci_list, args.out_dir, 
                      args.sample, args.flank, args.max_hits)