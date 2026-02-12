import os
import argparse
from Bio import SeqIO

def gather_sequences(base_dir, output_file, samples, type_dir, gene_id, min_taxa=4):
    """
    Gathers scattered gene files and combines them into one.
    """
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    combined_records = []
    found_samples = []

    # 1. Traverse each sample folder and find the corresponding gene file
    for sample in samples:
        # Path: results/05_extracted/{Sample}/{Type}/{Gene}.fasta
        gene_path = os.path.join(base_dir, sample, type_dir, f"{gene_id}.fasta")
        
        if os.path.exists(gene_path):
            try:
                # Read file
                records = list(SeqIO.parse(gene_path, "fasta"))
                if records:
                    combined_records.extend(records)
                    found_samples.append(sample)
            except Exception as e:
                print(f"⚠️ Error reading {gene_path}: {e}")

    # 2. Filtering: Are there at least N minimum samples gathered?
    unique_taxa_count = len(set(found_samples))
    
    if unique_taxa_count >= min_taxa:
        # Save file
        with open(output_file, "w") as f:
            SeqIO.write(combined_records, f, "fasta")
        print(f"✅ Gathered {gene_id} ({type_dir}): {unique_taxa_count} samples found.")
    else:
        # Create empty file if not enough samples gathered (error prevention)
        print(f"⚠️ Skipped {gene_id} ({type_dir}): Only {unique_taxa_count} samples (min={min_taxa}).")
        open(output_file, 'w').close()
        
        # Record skipped gene ID to a separate file
        skipped_log_file = os.path.join(os.path.dirname(output_file), "skipped_genes.txt")
        with open(skipped_log_file, "a") as log:
            log.write(f"{gene_id}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--samples", nargs='+', required=True)
    parser.add_argument("--type_dir", required=True) # Best or All
    parser.add_argument("--gene_id", required=True)
    parser.add_argument("--min_taxa", type=int, default=4)
    
    args = parser.parse_args()
    
    gather_sequences(
        base_dir=args.base_dir,
        output_file=args.out_file,
        samples=args.samples,
        type_dir=args.type_dir,
        gene_id=args.gene_id,
        min_taxa=args.min_taxa
    )