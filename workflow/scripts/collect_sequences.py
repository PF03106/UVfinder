import os
import argparse
from Bio import SeqIO

def collect_sequences(loci_list_path, input_dirs, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    with open(loci_list_path, 'r') as f:
        target_loci = [line.strip() for line in f if line.strip()]

    for locus in target_loci:
        all_records = []
        for s_dir in input_dirs:
            file_path = os.path.join(s_dir, f"{locus}.fasta")
            if os.path.exists(file_path):
                # QC 없이 모든 레코드를 그대로 추가합니다.
                all_records.extend(list(SeqIO.parse(file_path, "fasta")))

        if all_records:
            out_path = os.path.join(output_dir, f"{locus}_all.fasta")
            SeqIO.write(all_records, out_path, "fasta")
            
    print(f"✅ Collection complete for {len(target_loci)} loci.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci_list", required=True)
    parser.add_argument("--input_dirs", nargs="+", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()
    collect_sequences(args.loci_list, args.input_dirs, args.output_dir)