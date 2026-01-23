#!/usr/bin/env python3
import os
import re
import sys

def clean_original_id(header_line):
    """
    Extracts Chromosome or Scaffold numbers from complex FASTA headers.
    If no pattern is found, it cleans the header string for bioinformatics safety.
    """
    # Remove '>' and leading/trailing whitespaces
    full_header = header_line.strip()[1:]
    
    # 1. Try to find Chromosome patterns (e.g., Chromosome 1, Chr_02, LG1)
    chrom_match = re.search(r'chromosome\s*:?\s*([0-9]+)', full_header, re.IGNORECASE) \
               or re.search(r'chromosome([0-9]+)', full_header, re.IGNORECASE) \
               or re.search(r'chr[ _-]?([0-9]+)', full_header, re.IGNORECASE) \
               or re.search(r'LG([A-Za-z0-9]+)', full_header, re.IGNORECASE)
    
    if chrom_match:
        return f"Chr{chrom_match.group(1)}"

    # 2. Try to find Scaffold patterns (e.g., scaffold_69, Scaffold10, scaffold 14)
    scaff_match = re.search(r'scaffold[ _-]?([0-9]+)', full_header, re.IGNORECASE) \
               or re.search(r'scaffold\s*([0-9]+)', full_header, re.IGNORECASE)
    
    if scaff_match:
        # Standardize to 'Scaffold' + Number
        return f"Scaffold{scaff_match.group(1)}"

    # 3. Fallback: If no patterns match, use the first word (often Accession) and clean it
    # Take only the first word before a space
    primary_id = full_header.split()[0]
    # Replace non-alphanumeric chars with underscore
    cleaned = re.sub(r'[^a-zA-Z0-9]', '_', primary_id)
    cleaned = re.sub(r'_+', '_', cleaned).strip('_')
    
    return cleaned

def rename_fasta_headers(input_fna, output_fna, mapping_fna, s_id, order, genus, species):
    """
    Main function to process FASTA files and generate mapping logs.
    """
    # Standardized prefix for the pipeline
    prefix = f"{s_id}_{order}_{genus}_{species}"
    
    with open(input_fna, "r", encoding="utf-8", errors="ignore") as f_in, \
         open(output_fna, "w", encoding="utf-8") as f_out, \
         open(mapping_fna, "w", encoding="utf-8") as f_map:
        
        f_map.write("original_id\tnew_id\n")
        
        for line in f_in:
            if line.startswith(">"):
                original_header = line.strip()
                # Use the enhanced cleaning logic
                new_clean_id = clean_original_id(original_header)
                new_header = f">{prefix}_{new_clean_id}"
                
                f_out.write(new_header + "\n")
                # Store full original header for mapping traceability
                f_map.write(f"{original_header[1:]}\t{new_header[1:]}\n")
            else:
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 rename_fasta_headers.py <in.fna> <out.fna> <map.tsv> <ID> <Order> <Genus> <Species>")
        sys.exit(1)

    input_path, output_path, mapping_path = sys.argv[1:4]
    s_id, order_name, genus_name, species_name = sys.argv[4:8]

    print(f"🚀 Processing: {s_id} - Using logic for Chrs and Scaffolds")
    rename_fasta_headers(input_path, output_path, mapping_path, s_id, order_name, genus_name, species_name)
    print(f"✅ Renaming complete. Log saved to {mapping_path}")