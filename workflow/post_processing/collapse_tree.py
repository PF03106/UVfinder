import glob
import os
import argparse
from ete3 import Tree

# 1. Setup argument parser to accept directory path from the command line
parser = argparse.ArgumentParser(description="Collapse low bootstrap support nodes in tree files.")
parser.add_argument("-d", "--dir", required=True, help="Target directory containing .treefile files")
parser.add_argument("-t", "--threshold", type=float, default=10, help="Support threshold for collapsing nodes (default: 10)")
args = parser.parse_args()

target_dir = args.dir
THRESHOLD = args.threshold

# 2. Retrieve all .treefile paths from the specified directory
search_pattern = os.path.join(target_dir, "*.treefile")
tree_files = glob.glob(search_pattern)

print(f"Processing {len(tree_files)} trees in [{target_dir}]...\n")

for tree_file in tree_files:
    t = Tree(tree_file)
    # 3. Traverse the tree and collapse nodes with support values below the threshold
    for node in t.get_descendants():
        # Collapse internal nodes that do not meet the bootstrap support threshold
        if not node.is_leaf() and node.support < THRESHOLD:
            node.delete()
            
    # 4. Define output path and filename
    dir_name = os.path.dirname(tree_file)
    base_name = os.path.basename(tree_file)
    output_name = os.path.join(dir_name, f"{int(THRESHOLD)}_collapsed_{base_name}")
    
    # 5. Save the processed tree
    t.write(outfile=output_name)
    print(f"Successfully processed: {base_name} -> {os.path.basename(output_name)}")