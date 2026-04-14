#!/usr/bin/env python3

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

    count_total = 0
    locus_matched = 0
    locus_fallback = 0
    fallback_loci_list=[]

    # 3. Traverse probe directory and filter locus by locus
    with open(output_fasta, "w") as out_f:
        if not os.path.exists(probe_dir):
            print(f"❌ Error: Probe directory not found: {probe_dir}")
            return

        for locus_file in os.listdir(probe_dir):
            if not locus_file.endswith(".fasta"): continue
            
            # Extract locus ID from filename (e.g.: G4471.fasta -> G4471)
            locus_id = locus_file.split(".")[0]
            locus_path = os.path.join(probe_dir, locus_file)
            
            # load all sequences from the current locus file into memory
            all_records = list(SeqIO.parse(locus_path, "fasta"))
            
            # Filter sequences that start with the target order (case-insensitive)
            matched_records = [r for r in all_records if r.id.lower().startswith(order_norm)]
            
            # If there are matched records, use them; otherwise, use all records to avoid missing probes for this locus
            if matched_records:
                records_to_write = matched_records
                locus_matched += 1
            else:
                records_to_write = all_records
                locus_fallback += 1
                fallback_loci_list.append(locus_id)

            # Write the selected records to the output FASTA file, modifying the header to include locus ID for easier parsing later
            for record in records_to_write:
                # [Important] Change to 'LocusID|OriginalHeader' format for easier parsing later
                record.id = f"{locus_id}|{record.id}"
                record.description = record.id 
                
                SeqIO.write(record, out_f, "fasta")
                count_total += 1

    # 4. Summary (log file)
    print(f"✅ Successfully collected {count_total} total probes for {sample_id}.")
    print(f"   - Loci with target order '{order_norm}': {locus_matched}")
    print(f"   - Loci using all sequences (fallback): {locus_fallback}")
    if fallback_loci_list:
        print(f"   - Fallback loci: {', '.join(fallback_loci_list)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--samples_tsv", required=True, help="Path to samples.tsv")
    parser.add_argument("--probe_dir", required=True, help="Directory containing probe FASTA files")
    parser.add_argument("--output", required=True, help="Output FASTA file path")
    
    args = parser.parse_args()
    
    filter_probes_by_order(args.sample, args.samples_tsv, args.probe_dir, args.output)