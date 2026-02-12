import sys
from Bio import SeqIO

def filter_alignment(input_aln, output_aln, sample_list):
    """
    Keep only sequences starting with IDs from the Sample List in the Alignment file
    (Purpose: Remove Reference sequences)
    """
    kept_sequences = []
    
    # 1. Read FASTA file
    # MAFFT output is in FASTA format
    for record in SeqIO.parse(input_aln, "fasta"):
        is_my_sample = False
        for s in sample_list:
            if record.id.startswith(s):
                is_my_sample = True
                break
        
        if is_my_sample:
            kept_sequences.append(record)

    # 2. Save
    if kept_sequences:
        with open(output_aln, "w") as out_handle:
            SeqIO.write(kept_sequences, out_handle, "fasta")
        print(f"✅ Filtered: Kept {len(kept_sequences)} sequences (Removed Reference).")
    else:
        print("⚠️ Warning: No sequences kept after filtering!")
        # Create an empty file
        open(output_aln, 'w').close()

if __name__ == "__main__":
    # python script.py input.aln output.aln S001 S002 ...
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    samples = sys.argv[3:]
    
    filter_alignment(input_file, output_file, samples)