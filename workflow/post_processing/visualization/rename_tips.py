import sys
import re

def main():
    if len(sys.argv) != 4:
        print("Usage: python rename_tips.py [tree_file] [metadata_file] [output_file]")
        sys.exit(1)

    tree_file = sys.argv[1]
    meta_file = sys.argv[2]
    out_file = sys.argv[3]

    # 1. Read metadata (TSV) and create a mapping dictionary
    id_to_name = {}
    try:
        with open(meta_file, 'r') as f:
            header = f.readline()  # Skip the header row
            for line in f:
                parts = line.strip().split('\t')
                # Process lines with at least 4 columns: id, order, genus, species
                if len(parts) >= 4:
                    s_id = parts[0].strip()
                    genus = parts[2].strip()
                    species = parts[3].strip()
                    # Combine into Genus_species format (e.g., Antitrichia_curtipendula)
                    new_name = f"{genus}_{species}"
                    id_to_name[s_id] = new_name
    except Exception as e:
        print(f"Error reading metadata file: {e}")
        sys.exit(1)

    # 2. Read the original tree file
    try:
        with open(tree_file, 'r') as f:
            tree_text = f.read()
    except Exception as e:
        print(f"Error reading tree file: {e}")
        sys.exit(1)

    # 3. Use Regular Expressions to safely replace IDs in both formats
    # Sort IDs by length in descending order to prevent partial matching (e.g., matching 'S001' inside 'S0011')
    for s_id in sorted(id_to_name.keys(), key=len, reverse=True):
        new_name = id_to_name[s_id]
        # Regex explanation:
        # \b{s_id}(?=\b|_) matches 's_id' only if it is followed by a word boundary (\b) or an underscore (_)
        # This handles both standalone cases (S0001) and concatenated cases (S0001_G5280...)
        tree_text = re.sub(rf'\b{s_id}(?=\b|_)', new_name, tree_text)

    # 4. Save the modified tree text to the output file
    try:
        with open(out_file, 'w') as f:
            f.write(tree_text)
        print(f"✅ Success! Total of {len(id_to_name)} IDs processed.")
        print(f"✅ Output file saved to: {out_file}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()