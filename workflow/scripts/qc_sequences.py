import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def check_translation_quality(dna_record):
    """
    Checks if the DNA sequence contains internal stop codons when translated.
    """
    dna_seq = dna_record.seq
    
    # We use the standard genetic code (table 1).
    # In a more advanced version, we could check all 3 frames based on BLAST coordinates.
    try:
        aa_seq = dna_seq.translate()
        # Internal stop codons (denoted by '*') indicate frameshifts or poor assembly.
        # We allow a stop codon at the very end.
        if "*" in str(aa_seq)[:-1]:
            return False
        return True
    except:
        return False

def collect_and_qc(loci_list_path, input_dirs, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # Load interesting loci IDs
    with open(loci_list_path, 'r') as f:
        target_loci = [line.strip() for line in f if line.strip()]

    for locus in target_loci:
        collected_records = []
        
        for s_dir in input_dirs:
            file_path = os.path.join(s_dir, f"{locus}.fasta")
            
            if os.path.exists(file_path):
                for record in SeqIO.parse(file_path, "fasta"):
                    # Quality Check: Remove sequences with internal stops
                    if check_translation_quality(record):
                        collected_records.append(record)
                    else:
                        print(f"Skipping {record.id} in {locus} due to internal stop codon.")

        # Save merged sequences for this locus if we have results from at least one sample
        if collected_records:
            out_path = os.path.join(output_dir, f"{locus}_merged.fasta")
            SeqIO.write(collected_records, out_path, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci_list", required=True)
    parser.add_argument("--input_dirs", nargs="+", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()
    
    collect_and_qc(args.loci_list, args.input_dirs, args.output_dir)