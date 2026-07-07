#!/usr/bin/env python3
"""
astral_triple_plot.py  –  Plot 3 wASTRAL trees aligned with maximized layout rotation and branch support
"""
import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import re
import sys
import argparse
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces

DEFAULT_COLORS = ["#4D94FF", "#FF4D4D", "#FFD700"]

def load_tree(tree_path: str) -> Tree:
    with open(tree_path) as fh:
        nwk = fh.read().strip()
    return Tree(nwk, format=1, quoted_node_names=True)

_Q_PATTERN = re.compile(r'q([123])=([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)')
_PP_PATTERN = re.compile(r'pp1=([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)')
_SUPPORT_PATTERN = re.compile(r'support=([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)')

def parse_annotation(name: str) -> dict:
    result = {}
    for m in _Q_PATTERN.finditer(name):
        result[f"q{m.group(1)}"] = float(m.group(2))
    
    m2 = _PP_PATTERN.search(name)
    if m2:
        result["pp1"] = float(m2.group(1))
        
    m3 = _SUPPORT_PATTERN.search(name)
    if m3:
        result["support"] = float(m3.group(1))
        
    return result

def root_on_outgroup(tree: Tree, outgroup_names: list) -> Tree:
    leaves = {leaf.name: leaf for leaf in tree.get_leaves()}
    for n in outgroup_names:
        if n not in leaves:
            raise ValueError(f"Outgroup {n} not found.")
    if len(outgroup_names) == 1:
        tree.set_outgroup(leaves[outgroup_names[0]])
    else:
        tree.set_outgroup(tree.get_common_ancestor([leaves[n] for n in outgroup_names]))
    return tree

def extract_split_annotations(tree: Tree) -> dict:
    split_map = {}
    for node in tree.traverse():
        if not node.is_leaf() and node.name:
            leaf_set = frozenset(leaf.name for leaf in node.get_leaves())
            split_map[leaf_set] = node.name
    return split_map

def reassign_annotations(tree: Tree, split_map: dict) -> None:
    all_leaves = frozenset(leaf.name for leaf in tree.get_leaves())
    for node in tree.traverse():
        if node.is_leaf():
            continue
        leaf_set = frozenset(leaf.name for leaf in node.get_leaves())
        ann = split_map.get(leaf_set) or split_map.get(all_leaves - leaf_set)
        node.name = ann if ann else ""

# 💡 Flip
def align_tree_to_target(tree: Tree, target_order: list):
    target_idx = {name: i for i, name in enumerate(target_order)}

    def get_avg_index(n):
        valid_indices = [target_idx[l.name] for l in n.get_leaves() if l.name in target_idx]
        return sum(valid_indices) / len(valid_indices) if valid_indices else float('inf')

    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.children.sort(key=get_avg_index)

def apply_styling(tree: Tree, colors: list, pie_size: int, show_pp: bool, show_support: bool, branch_lw: float, tip_font: int):
    for node in tree.traverse():
        ns = NodeStyle()
        ns["size"] = 0
        ns["hz_line_width"] = branch_lw
        ns["vt_line_width"] = branch_lw
        node.set_style(ns)

        if node.is_leaf():
            tf = TextFace(" " + node.name, fsize=tip_font, fgcolor="#171717")
            tf.margin_left = 4
            node.add_face(tf, column=0, position="branch-right")
        else:
            ann = parse_annotation(node.name) if node.name else {}
            if ann:
                q1, q2, q3 = ann.get("q1", 0), ann.get("q2", 0), ann.get("q3", 0)
                total = q1 + q2 + q3
                if total > 0:
                    percents = [100 * q1 / total, 100 * q2 / total, 100 * q3 / total]
                    pie = faces.PieChartFace(percents, colors=colors, width=pie_size, height=pie_size)
                    pie.border.width = 0.5
                    node.add_face(pie, column=0, position="branch-top")
                
                if show_pp and "pp1" in ann:
                    pp_face = TextFace(f" {ann['pp1']:.2f}", fsize=7, fgcolor="#171717")
                    node.add_face(pp_face, column=1, position="branch-top")
                    
                if show_support and "support" in ann:
                    sup_face = TextFace(f" {ann['support']:.2f}", fsize=6, bold=True, fgcolor="#C00000")
                    node.add_face(sup_face, column=1, position="branch-bottom")

def build_legend(ts: TreeStyle, colors: list):
    labels = ["q1 (focal)", "q2 (alt. 1)", "q3 (alt. 2)"]
    for color, label in zip(colors, labels):
        ts.legend.add_face(faces.RectFace(12, 12, color, color), column=0)
        ts.legend.add_face(TextFace(f"  {label}", fsize=9), column=1)
    ts.legend_position = 4

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-t", "--trees", required=True, help="3 Newick tree files, comma-separated")
    p.add_argument("-o", "--output", default="astral_triple_aligned.svg")
    p.add_argument("-g", "--outgroup", default=None, help="Outgroup taxon name")
    p.add_argument("--pie-size", type=int, default=20)
    p.add_argument("--scale", type=float, default=50)
    p.add_argument("--show-support", action="store_true", help="Display branch support values")
    args = p.parse_args()

    tree_files = [f.strip() for f in args.trees.split(",")]
    if len(tree_files) != 3:
        sys.exit("Error: Please provide exactly 3 tree files separated by commas.")

    # 1. Load, Root, and Restore
    trees = []
    for f in tree_files:
        t = load_tree(f)
        split_map = extract_split_annotations(t)
        if args.outgroup:
            og_names = [n.strip() for n in args.outgroup.split(",")]
            t = root_on_outgroup(t, og_names)
            reassign_annotations(t, split_map)
        trees.append(t)

    # 2. 첫 번째 트리 기준으로 나머지 트리들 강제 회전(Flipping) 정렬
    trees[0].ladderize()
    target_order = [leaf.name for leaf in trees[0].get_leaves()]

    for t in trees[1:]:
        t.ladderize()
        align_tree_to_target(t, target_order)

    # 3. Set Style
    for i, t in enumerate(trees):
        apply_styling(t, DEFAULT_COLORS, args.pie_size, True, args.show_support, 2.0, 10)

    # 4. Set canvas
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.scale = args.scale
    ts.branch_vertical_margin = 10
    ts.force_topology = False
    build_legend(ts, DEFAULT_COLORS)

    ts.title.add_face(TextFace(" 3 Aligned Trees Comparison (q1/q2/q3 & Support)", fsize=10, bold=True), column=0)

    master = Tree()
    for t in trees:
        master.add_child(t)
        
    master.render(args.output, tree_style=ts, w=1400, units="px", dpi=150)
    print(f"🎉 Visulaized three trees: {args.output}")

if __name__ == "__main__":
    main()