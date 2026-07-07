import re
import os
import csv
import argparse
from ete3 import Tree

# Ensure the script runs in headless mode for HPC environments
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

def load_and_clean_tree(file_path):
    """Reads a Newick file and pre-processes metadata for ete3 compatibility."""
    try:
        if not os.path.exists(file_path):
            print(f"Error: File not found at {file_path}")
            return None
            
        with open(file_path, 'r') as f:
            raw_newick = f.read().replace("'", "")
            
            # Replace semicolons within brackets with pipes to prevent parsing errors
            def clean_meta(match):
                safe_inner = match.group(1).replace(";", "|")
                return f"META_{safe_inner}_END"
            
            cleaned_newick = re.sub(r'\[(.*?)\]', clean_meta, raw_newick)
            return Tree(cleaned_newick, format=1)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def get_quartet_map(tree):
    """Traverses the tree and maps clade members to their Quartet Scores."""
    score_map = {}
    for node in tree.traverse():
        if not node.is_leaf() and node.name:
            clade_key = frozenset(node.get_leaf_names())
            
            q1 = re.search(r'q1=([0-9.]+)', node.name)
            q2 = re.search(r'q2=([0-9.]+)', node.name)
            q3 = re.search(r'q3=([0-9.]+)', node.name)
            
            if q1 and q2 and q3:
                score_map[clade_key] = (
                    float(q1.group(1)), 
                    float(q2.group(1)), 
                    float(q3.group(1))
                )
    return score_map

def main():
    parser = argparse.ArgumentParser(description="Compare Quartet Scores between two ASTER trees and save to TSV.")
    parser.add_argument("-i1", "--input1", required=True, help="Path to the first tree file (.tre)")
    parser.add_argument("-i2", "--input2", required=True, help="Path to the second tree file (.tre)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory to save the TSV file")
    parser.add_argument("-n", "--name", default="quartet_comparison.tsv", help="Output TSV filename (default: quartet_comparison.tsv)")
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        print(f"Created output directory: {args.outdir}")

    print("Loading and processing trees...")
    tree1 = load_and_clean_tree(args.input1)
    tree2 = load_and_clean_tree(args.input2)

    if not tree1 or not tree2:
        print("Error: Tree loading failed. Exiting.")
        return

    map1 = get_quartet_map(tree1)
    map2 = get_quartet_map(tree2)

    output_path = os.path.join(args.outdir, args.name)

    common_clades = set(map1.keys()) & set(map2.keys())
    unique_to_1 = len(map1) - len(common_clades)
    unique_to_2 = len(map2) - len(common_clades)

    # 3. Write results to TSV
    try:
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write Header
            writer.writerow(["Clade_Members", "T1_q1", "T1_q2", "T1_q3", "T2_q1", "T2_q2", "T2_q3"])
            
            # Write Common Clades Data
            for clade in sorted(list(common_clades), key=len):
                leaf_names = ",".join(sorted(list(clade)))
                s1 = map1[clade]
                s2 = map2[clade]
                writer.writerow([leaf_names, s1[0], s1[1], s1[2], s2[0], s2[1], s2[2]])
            
            # Add Summary Statistics at the end
            writer.writerow([])
            writer.writerow(["Tree 1: " + os.path.basename(args.input1)])
            writer.writerow(["Tree 2: " + os.path.basename(args.input2)])
            writer.writerow(["SUMMARY_STATISTICS"])
            writer.writerow(["Common_Nodes", len(common_clades)])
            writer.writerow(["Unique_to_Tree1", unique_to_1])
            writer.writerow(["Unique_to_Tree2", unique_to_2])
            writer.writerow(["Total_Nodes_Tree1", len(map1)])
            writer.writerow(["Total_Nodes_Tree2", len(map2)])

        print(f"\nSuccessfully saved comparison to: {output_path}")
        print(f" - Common nodes found: {len(common_clades)}")
        print(f" - Unique to Tree 1: {unique_to_1}")
        print(f" - Unique to Tree 2: {unique_to_2}")

    except Exception as e:
        print(f"Error writing to TSV: {e}")

if __name__ == "__main__":
    main()