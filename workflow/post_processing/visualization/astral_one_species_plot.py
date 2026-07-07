#!/usr/bin/env python3
"""
astral_pie_plot.py  –  Plot weighted ASTRAL trees with q-value pie charts

Each internal node is annotated with a pie chart showing the three quartet
support values (q1, q2, q3) that weighted ASTRAL embeds in the node-label
field.  The pie chart convention follows ASTRAL-III / wASTRAL output:

    q1 = support for the displayed topology
    q2 = support for the first alternative
    q3 = support for the second alternative

Usage
-----
    python astral_pie_plot.py -t tree.nwk -g Yarrowia_lipolytica -o out.svg
    python astral_pie_plot.py -t tree.nwk -g "Taxon_A,Taxon_B" -o out.svg
    python astral_pie_plot.py -n "((...));" -g Outgroup -o out.svg

Options
-------
    -t / --tree FILE          Newick tree file
    -n / --newick STRING      Newick string literal (alternative to -t)
    -o / --output FILE        Output file  [default: astral_tree.svg]
                              Extension controls format: .svg, .pdf, .png
    -g / --outgroup TAXA      Outgroup taxon name(s), comma-separated
    --pie-size INT            Pie chart diameter in pixels   [default: 22]
    --scale FLOAT             Pixels per branch-length unit  [default: 60]
    --width INT               Rendered width in pixels       [default: 900]
    --colors C1,C2,C3         Hex colours for q1/q2/q3       [default: blue/red/yellow]
    --show-lengths            Label branch lengths on the tree
    --show-pp                 Show pp1 (posterior prob.) as a number below each pie
    --branch-lw FLOAT         Branch line width               [default: 2]
    --tip-font INT            Tip label font size             [default: 11]
    --no-legend               Suppress the colour legend
    --midpoint                Root at midpoint (ignores -g)
"""

import re
import sys
import argparse
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces


# ── Colour palette ─────────────────────────────────────────────────────────────
DEFAULT_COLORS = ["#4D94FF", "#FF4D4D", "#FFD700"]   # q1 blue / q2 red / q3 yellow


# ── Newick parsing ──────────────────────────────────────────────────────────────
def load_tree(tree_arg: str, is_string: bool = False) -> Tree:
    """
    Load an ete3 Tree from a file path or a raw Newick string.
    format=1  →  internal node names are preserved (required for ASTRAL labels).
    """
    if is_string:
        return Tree(tree_arg, format=1, quoted_node_names=True)
    with open(tree_arg) as fh:
        nwk = fh.read().strip()
    return Tree(nwk, format=1, quoted_node_names=True)


# ── Q-value extraction ──────────────────────────────────────────────────────────
_Q_PATTERN = re.compile(r'q([123])=([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)')
_PP_PATTERN = re.compile(r'pp1=([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)')

def parse_annotation(name: str) -> dict:
    """Return a dict with q1, q2, q3 (and pp1 if present) from an ASTRAL label."""
    result = {}
    for m in _Q_PATTERN.finditer(name):
        result[f"q{m.group(1)}"] = float(m.group(2))
    m2 = _PP_PATTERN.search(name)
    if m2:
        result["pp1"] = float(m2.group(1))
    return result


# ── Rooting ─────────────────────────────────────────────────────────────────────
def root_on_outgroup(tree: Tree, outgroup_names: list) -> Tree:
    """
    Root the tree on an outgroup.  One name → set_outgroup on that leaf.
    Multiple names → root on the MRCA of those leaves.
    """
    leaves = {leaf.name: leaf for leaf in tree.get_leaves()}
    missing = [n for n in outgroup_names if n not in leaves]
    if missing:
        raise ValueError(f"Outgroup taxon/taxa not found in tree: {missing}\n"
                         f"Available names: {sorted(leaves)}")

    if len(outgroup_names) == 1:
        tree.set_outgroup(leaves[outgroup_names[0]])
    else:
        og_nodes = [leaves[n] for n in outgroup_names]
        mrca = tree.get_common_ancestor(og_nodes)
        tree.set_outgroup(mrca)
    return tree


# ── Annotation invariance across rerooting ─────────────────────────────────────
def extract_split_annotations(tree: Tree) -> dict:
    """
    Snapshot { frozenset(leaf_names) -> annotation_string } for every
    annotated internal node BEFORE rerooting.

    ASTRAL q-values describe edge bipartitions in the *unrooted* tree, which
    are fully defined by which leaf set is on each side of the edge.  That
    mapping is invariant under rerooting, so we can restore correct annotations
    after ete3 restructures the node objects.
    """
    split_map = {}
    for node in tree.traverse():
        if not node.is_leaf() and node.name:
            leaf_set = frozenset(leaf.name for leaf in node.get_leaves())
            split_map[leaf_set] = node.name
    return split_map


def reassign_annotations(tree: Tree, split_map: dict) -> None:
    """
    Restore ASTRAL annotations after rerooting.

    For each internal node, the correct annotation is the one whose original
    bipartition matches this node's current leaf set — or its complement,
    since a split {A|B} is the same bipartition viewed from either side.
    The tree root gets no annotation (no well-defined split above it).
    """
    all_leaves = frozenset(leaf.name for leaf in tree.get_leaves())
    for node in tree.traverse():
        if node.is_leaf():
            continue
        leaf_set = frozenset(leaf.name for leaf in node.get_leaves())
        ann = split_map.get(leaf_set) or split_map.get(all_leaves - leaf_set)
        node.name = ann if ann else ""


# ── Node styling ────────────────────────────────────────────────────────────────
def style_node(node, ann: dict, colors: list, pie_size: int,
               show_pp: bool, branch_lw: float):
    """Apply NodeStyle + optional PieChartFace to one internal node."""
    ns = NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = branch_lw
    ns["vt_line_width"] = branch_lw
    node.set_style(ns)

    if not ann:
        return

    q1, q2, q3 = ann.get("q1", 0), ann.get("q2", 0), ann.get("q3", 0)
    total = q1 + q2 + q3
    if total == 0:
        return

    percents = [100 * q1 / total, 100 * q2 / total, 100 * q3 / total]
    pie = faces.PieChartFace(
        percents,
        colors=colors,
        width=pie_size,
        height=pie_size,
    )
    pie.border.width = 0.5
    node.add_face(pie, column=0, position="float")

    if show_pp and "pp1" in ann:
        pp_face = TextFace(f" {ann['pp1']:.2f}", fsize=7, fgcolor="#333333")
        node.add_face(pp_face, column=1, position="branch-bottom")


def style_leaf(node, tip_font: int, branch_lw: float):
    ns = NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = branch_lw
    ns["vt_line_width"] = branch_lw
    node.set_style(ns)
    # ete3 renders leaf names automatically; optionally override size here
    # (TreeStyle.aligned_foot_header controls tip label font globally in
    #  older ete3 versions — we can also patch it per node):
    tf = TextFace(" " + node.name, fsize=tip_font, fgcolor="#111111")
    tf.margin_left = 4
    node.add_face(tf, column=0, position="branch-right")


# ── Legend ──────────────────────────────────────────────────────────────────────
def build_legend(ts: TreeStyle, colors: list):
    labels = ["q1 (focal topology)", "q2 (alt. 1)", "q3 (alt. 2)"]
    for color, label in zip(colors, labels):
        ts.legend.add_face(faces.RectFace(12, 12, color, color), column=0)
        ts.legend.add_face(TextFace(f"  {label}", fsize=9), column=1)
    ts.legend_position = 4   # bottom-right


# ── Main plotting logic ─────────────────────────────────────────────────────────
def plot(args):
    # 1. Load tree
    if args.newick:
        tree = load_tree(args.newick, is_string=True)
    elif args.tree:
        tree = load_tree(args.tree, is_string=False)
    else:
        sys.exit("Error: provide --tree or --newick.")

    # 2. Root
    # Snapshot bipartition->annotation BEFORE rerooting; ete3 moves node
    # objects around during set_outgroup(), scrambling which annotation sits
    # at which position.  We restore correct assignments by leaf-set lookup.
    split_map = extract_split_annotations(tree)

    if args.midpoint:
        tree.set_outgroup(tree.get_midpoint_outgroup())
        reassign_annotations(tree, split_map)
    elif args.outgroup:
        og_names = [n.strip() for n in args.outgroup.split(",")]
        tree = root_on_outgroup(tree, og_names)
        reassign_annotations(tree, split_map)

    # 3. Colours
    if args.colors:
        colors = [c.strip() for c in args.colors.split(",")]
        if len(colors) != 3:
            sys.exit("--colors requires exactly 3 comma-separated colour values.")
    else:
        colors = DEFAULT_COLORS

    pie_size  = args.pie_size
    branch_lw = args.branch_lw
    tip_font  = args.tip_font
    show_pp   = args.show_pp

    # 4. Apply styles
    for node in tree.traverse():
        if node.is_leaf():
            style_leaf(node, tip_font, branch_lw)
        else:
            ann = parse_annotation(node.name) if node.name else {}
            style_node(node, ann, colors, pie_size, show_pp, branch_lw)

    # 5. Global tree style
    ts = TreeStyle()
    ts.show_leaf_name   = False   # we add TextFace manually for better control
    ts.show_branch_length = args.show_lengths
    ts.show_branch_support = False
    ts.scale = args.scale
    ts.mode  = "r"                # rectangular
    ts.branch_vertical_margin = 8

    # Title
    ts.title.add_face(TextFace("  wASTRAL quartet support (q1/q2/q3)",
                                fsize=11, bold=True, fgcolor="#444444"), column=0)

    if not args.no_legend:
        build_legend(ts, colors)

    # 6. Render
    out = args.output
    tree.render(out, tree_style=ts, w=args.width, units="px", dpi=150)
    print(f"Tree written to: {out}")


# ── CLI ─────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Plot a weighted ASTRAL tree with q-value pie charts at nodes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument("-t", "--tree",   metavar="FILE",   help="Newick tree file")
    src.add_argument("-n", "--newick", metavar="STRING", help="Newick string literal")

    p.add_argument("-o", "--output",    default="astral_tree.svg", metavar="FILE")
    p.add_argument("-g", "--outgroup",  default=None,  metavar="TAXA",
                   help="Outgroup taxon name(s), comma-separated")
    p.add_argument("--midpoint",        action="store_true",
                   help="Root at midpoint (overrides --outgroup)")
    p.add_argument("--pie-size",  type=int,   default=22,   metavar="INT")
    p.add_argument("--scale",     type=float, default=60,   metavar="FLOAT")
    p.add_argument("--width",     type=int,   default=900,  metavar="INT")
    p.add_argument("--colors",    default=None, metavar="C1,C2,C3")
    p.add_argument("--show-lengths",  action="store_true", dest="show_lengths")
    p.add_argument("--show-pp",       action="store_true", dest="show_pp")
    p.add_argument("--branch-lw", type=float, default=2.0,  dest="branch_lw")
    p.add_argument("--tip-font",  type=int,   default=11,   dest="tip_font")
    p.add_argument("--no-legend", action="store_true",      dest="no_legend")

    args = p.parse_args()
    plot(args)


if __name__ == "__main__":
    main()
