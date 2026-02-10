import os
import argparse
import pandas as pd
from Bio import SeqIO

def filter_probes_by_order(sample_id, samples_tsv, probe_dir, output_fasta):
    # 1. Load sample info
    try:
        df = pd.read_csv(samples_tsv, sep='\t').set_index("sample_id")
        # Get Order info for this sample
        sample_order = df.loc[sample_id, "order"]
    except KeyError:
        print(f"❌ Error: Sample ID '{sample_id}' not found in {samples_tsv}.")
        raise
    except Exception as e:
        print(f"❌ Error reading samples.tsv: {e}")
        raise

    # 2. Normalize order string (lowercase, remove spaces)
    # e.g.: "Bryales " -> "bryales"
    order_norm = str(sample_order).lower().replace(" ", "").strip()
    print(f"--- 🧬 Extracting probes for {sample_id} (Target Order: {order_norm}) ---")

    count = 0
    # 3. Traverse probe directory and filter
    with open(output_fasta, "w") as out_f:
        if not os.path.exists(probe_dir):
            print(f"❌ Error: Probe directory not found: {probe_dir}")
            return

        for locus_file in os.listdir(probe_dir):
            if not locus_file.endswith(".fasta"): continue
            
            # Extract locus ID from filename (e.g.: G4471.fasta -> G4471)
            locus_id = locus_file.split(".")[0]
            locus_path = os.path.join(probe_dir, locus_file)
            
            # Parse FASTA
            for record in SeqIO.parse(locus_path, "fasta"):
                # Check if header starts with the target order (e.g.: >Bryales_...)
                if record.id.lower().startswith(order_norm):
                    # [Important] Change to 'LocusID|OriginalHeader' format for easier parsing later
                    record.id = f"{locus_id}|{record.id}"
                    record.description = record.id 
                    
                    SeqIO.write(record, out_f, "fasta")
                    count += 1
    
    if count == 0:
        print(f"⚠️ Warning: No matching probes found for order '{order_norm}' in {sample_id}.")
    else:
        print(f"✅ Successfully collected {count} probes for {sample_id}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--samples_tsv", required=True, help="Path to samples.tsv")
    parser.add_argument("--probe_dir", required=True, help="Directory containing probe FASTA files")
    parser.add_argument("--output", required=True, help="Output FASTA file path")
    
    args = parser.parse_args()
    
    filter_probes_by_order(args.sample, args.samples_tsv, args.probe_dir, args.output)