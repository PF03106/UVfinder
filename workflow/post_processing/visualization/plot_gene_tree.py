#!/usr/bin/env python3
import os
# If you want to render on a headless server without GUI, keep this line active.
os.environ['QT_QPA_PLATFORM'] = 'offscreen' 

import argparse
import sys
import csv
from ete3 import Tree, TreeStyle, NodeStyle

def main():
    # 1. Parse command-line arguments
    parser = argparse.ArgumentParser(description="Root a tree using multiple outgroups, rename leaves using a TSV mapping, and render as SVG.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input Newick file (.tre)")
    # Use nargs='+' to accept one or multiple outgroup IDs as a list (e.g., -g S0003 S0005)
    parser.add_argument("-g", "--outgroup", required=True, nargs='+', help="Original sample IDs of the outgroups (e.g., S0003 S0005)")
    parser.add_argument("-m", "--meta", default="config/samples.tsv", help="Path to the metadata TSV file (e.g., config/samples.tsv)")
    parser.add_argument("-o", "--output", default="rooted_tree.svg", help="Output image file path (default: rooted_tree.svg)")
    args = parser.parse_args()

    # 2. Load metadata and build a mapping dictionary (ID -> Genus_species)
    id_to_name = {}
    try:
        with open(args.meta, 'r', encoding='utf-8') as f:
            # Use csv.DictReader to automatically map columns by their header names
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Combine 'genus' and 'species' columns with an underscore
                new_name = f"{row['genus']}_{row['species']}"
                id_to_name[row['sample_id']] = new_name
        print(f"✔ Successfully loaded mapping for {len(id_to_name)} samples.")
    except Exception as e:
        print(f"❌ Error reading metadata file: {e}")
        sys.exit(1)

    # 3. Load the tree from the Newick file
    try:
        t = Tree(args.input)
    except Exception as e:
        print(f"❌ Error reading the tree file: {e}")
        sys.exit(1)

    # 4. Root the tree using the target outgroup nodes (Must be done before renaming)
    target_nodes = []
    for name in args.outgroup:
        # Search for the node by its exact name in the tree
        found = t.search_nodes(name=name)
        if found:
            target_nodes.append(found[0])
        else:
            print(f"⚠ Warning: Outgroup node '{name}' not found in the tree.")

    if not target_nodes:
        print("❌ Error: None of the specified outgroup nodes were found.")
        sys.exit(1)

    try:
        if len(target_nodes) == 1:
            # If only one outgroup is provided, root it directly
            t.set_outgroup(target_nodes[0])
        else:
            # If multiple outgroups are provided, find their Lowest Common Ancestor (LCA)
            lca = t.get_common_ancestor(target_nodes)
            t.set_outgroup(lca)
        
        rooted_names = ', '.join([n.name for n in target_nodes])
        print(f"✔ Successfully rooted tree with outgroup(s): {rooted_names}")
    except Exception as e:
        print(f"❌ Error during rooting: {e}")
        sys.exit(1)

    # 5. Rename leaves based on the mapping dictionary
    for leaf in t.get_leaves():
        # Split the leaf name by underscore
        parts = leaf.name.split("_")
        leaf_id = parts[0]
        
        # If the sample ID exists in the metadata mapping, replace it
        if leaf_id in id_to_name:
            # Preserve the original suffix (e.g., rank, sex) if it exists
            if len(parts) > 1:
                suffix = "_".join(parts[1:])
                leaf.name = f"{id_to_name[leaf_id]}_{suffix}"
            else:
                leaf.name = id_to_name[leaf_id]

    # 6. Configure tree visualization styles (TreeStyle & NodeStyle)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = True # Display branch support values
    
    # Adjust layout settings for better readability
    ts.branch_vertical_margin = 15 
    
    # Customize node style for a clean, academic look
    nstyle = NodeStyle()
    nstyle["size"] = 0 # Remove default node circles
    nstyle["vt_line_width"] = 2 # Vertical branch line width
    nstyle["hz_line_width"] = 2 # Horizontal branch line width
    nstyle["vt_line_color"] = "#333333" # Dark gray lines
    nstyle["hz_line_color"] = "#333333"

    # Apply the node style to all nodes in the tree
    for n in t.traverse():
        n.set_style(nstyle)

    # 7. Render and save the tree image as SVG
    try:
        # Render the tree to a file (width set to 183mm)
        t.render(args.output, w=183, units="mm", tree_style=ts)
        print(f"✔ Tree rendered and saved to: {args.output}")
    except Exception as e:
        print(f"❌ Error rendering tree: {e}")

if __name__ == "__main__":
    main()