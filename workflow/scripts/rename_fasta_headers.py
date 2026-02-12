#!/usr/bin/env python3
import os
import re
import sys

def clean_original_id(header_line):
    """
    1. Standardizes Chromosome/LG headers (e.g., Chr1, ChrX).
    2. For other sequences (Scaffolds, etc.), preserves the full header but replaces 
       special characters with '_' for safety/compatibility.
    """
    # Remove '>' and leading/trailing white spaces
    full_header = header_line.strip()[1:]
    
    # --- [1] Chromosome ---
    chrom_match = re.search(r'(?:chromosome|chr|LG|gi)\s*[-_:|]?\s*([0-9]+|[UVuv]|Un)\b', full_header, re.IGNORECASE)
    
    if chrom_match:
        return f"Chr{chrom_match.group(1)}"

    # --- [2] else(Scaffold, etc.) ---

    cleaned = re.sub(r'[ \t,;|]+', '_', full_header)
    cleaned = re.sub(r'[^a-zA-Z0-9._-]', '', cleaned)
    cleaned = re.sub(r'_+', '_', cleaned).strip('_')
    
    return cleaned

def rename_fasta_headers(input_fna, output_fna, mapping_fna, s_id, order, genus, species):
    # prefix based on sample.tsv e.g.: S0028_dicriales_ceratodon_purpureus
    prefix = f"{s_id}_{order}_{genus}_{species}"
    
    with open(input_fna, "r", encoding="utf-8", errors="ignore") as f_in, \
         open(output_fna, "w", encoding="utf-8") as f_out, \
         open(mapping_fna, "w", encoding="utf-8") as f_map:
        
        f_map.write("original_id\tnew_id\n")
        
        for line in f_in:
            if line.startswith(">"):
                original_header = line.strip()
                new_clean_id = clean_original_id(original_header)
                new_header = f">{prefix}_{new_clean_id}"
                
                f_out.write(new_header + "\n")
                # Save mapping log
                f_map.write(f"{original_header[1:]}\t{new_header[1:]}\n")
            else:
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 rename_fasta_headers.py <in.fna> <out.fna> <map.tsv> <ID> <Order> <Genus> <Species>")
        sys.exit(1)

    input_path, output_path, mapping_path = sys.argv[1:4]
    s_id, order_name, genus_name, species_name = sys.argv[4:8]

    print(f"🚀 Processing: {s_id} - Safe renaming mode (Keep full headers)")
    rename_fasta_headers(input_path, output_path, mapping_path, s_id, order_name, genus_name, species_name)
    print(f"✅ Renaming complete. Log saved to {mapping_path}")