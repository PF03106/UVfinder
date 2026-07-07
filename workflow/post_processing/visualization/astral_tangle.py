#!/usr/bin/env python3
"""
astral_dual_face.py  –  Plot 2 wASTRAL trees facing each other (Tanglegram style)
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
    
    # 💡 파싱 에러 방지용 특수문자 청소 로직 (완벽하게 유지됨)
    nwk = nwk.replace("'", "")
    nwk = re.sub(r'\[(.*?)\]', lambda m: m.group(0).replace(';', '|'), nwk)
    return Tree(nwk, format=1)

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

# 💡 기준 트리의 팁 순서에 맞춰 나뭇가지를 회전(Flip)시키는 커스텀 함수
def align_tree_to_target(tree: Tree, target_order: list):
    target_idx = {name: i for i, name in enumerate(target_order)}
    
    def get_avg_index(n):
        valid_indices = [target_idx[l.name] for l in n.get_leaves() if l.name in target_idx]
        return sum(valid_indices) / len(valid_indices) if valid_indices else float('inf')

    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.children.sort(key=get_avg_index)

# 💡 [추가] is_right_tree 파라미터를 추가하여 텍스트 라벨 위치를 좌/우로 조절합니다.
def apply_styling(tree: Tree, colors: list, pie_size: int, show_pp: bool, show_support: bool, branch_lw: float, tip_font: int, is_right_tree: bool = False):
    for node in tree.traverse():
        ns = NodeStyle()
        ns["size"] = 0
        ns["hz_line_width"] = 1
        ns["vt_line_width"] = 1
        node.set_style(ns)

        if node.is_leaf():
            # 오른쪽 나무는 글자를 팁의 '왼쪽'으로 보내기 위해 여백을 줍니다.
            if is_right_tree:
                tf = TextFace(node.name + "   ", fsize=tip_font, fgcolor="#111111")
            else:
                tf = TextFace("   " + node.name, fsize=tip_font, fgcolor="#111111")
            node.add_face(tf, column=0, position="branch-right")
        else:
            ann = parse_annotation(node.name) if node.name else {}
            if ann:
                q1, q2, q3 = ann.get("q1", 0), ann.get("q2", 0), ann.get("q3", 0)
                total = q1 + q2 + q3
                if total > 0:
                    percents = [100 * q1 / total, 100 * q2 / total, 100 * q3 / total]
                    pie = faces.PieChartFace(percents, colors=colors, width=pie_size, height=pie_size)
                    pie.border.width = 0
                    node.add_face(pie, column=0, position="float")
                
                if show_pp and "pp1" in ann:
                    pp_face = TextFace(f" {ann['pp1']:.2f}", fsize=7, fgcolor="#333333")
                    node.add_face(pp_face, column=1, position="branch-bottom")
                    
                if show_support and "support" in ann:
                    sup_face = TextFace(f" {ann['support']:.2f}", fsize=9, bold=True, fgcolor="#000000")
                    node.add_face(sup_face, column=1, position="branch-top")

def build_legend(ts: TreeStyle, colors: list):
    labels = ["q1 (focal)", "q2 (alt. 1)", "q3 (alt. 2)"]
    for color, label in zip(colors, labels):
        ts.legend.add_face(faces.RectFace(12, 12, color, color), column=0)
        ts.legend.add_face(TextFace(f"  {label}", fsize=9), column=1)
    ts.legend_position = 4

def main():
    p = argparse.ArgumentParser()
    # 💡 트리를 딱 2개만 받도록 수정
    p.add_argument("-t", "--trees", required=True, help="2 Newick tree files, comma-separated (e.g., left.tre,right.tre)")
    p.add_argument("-o", "--output", default="astral_dual_face.svg")
    p.add_argument("-g", "--outgroup", default=None, help="Outgroup taxon name")
    p.add_argument("--pie-size", type=int, default=20)
    p.add_argument("--scale", type=float, default=50)
    p.add_argument("--show-support", action="store_true", help="Display branch support values")
    args = p.parse_args()

    tree_files = [f.strip() for f in args.trees.split(",")]
    if len(tree_files) != 2:
        sys.exit("Error: Please provide exactly 2 tree files separated by commas.")

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

    # 2. 첫 번째 트리(왼쪽) 기준으로 두 번째 트리(오른쪽) 강제 회전 정렬
    trees[0].ladderize()
    target_order = [leaf.name for leaf in trees[0].get_leaves()]

    trees[1].ladderize()
    align_tree_to_target(trees[1], target_order)

    # 3. 스타일 입히기 (오른쪽 트리는 is_right_tree=True 로 설정)
    apply_styling(trees[0], DEFAULT_COLORS, args.pie_size, False, args.show_support, 2.0, 10, is_right_tree=False)
    apply_styling(trees[1], DEFAULT_COLORS, args.pie_size, False, args.show_support, 2.0, 10, is_right_tree=True)

    # 4. [핵심] 각각의 트리를 독립적인 그림(TreeFace)으로 렌더링 설정
    # --- 왼쪽 트리 설정 (기본 방향) ---
    ts_left = TreeStyle()
    ts_left.show_leaf_name = False
    ts_left.show_branch_support = False
    ts_left.scale = args.scale
    ts_left.branch_vertical_margin = 10
    ts_left.orientation = 0 # 0 = 왼쪽에서 오른쪽으로 자람
    build_legend(ts_left, DEFAULT_COLORS)
    ts_left.title.add_face(TextFace(" Tree 1 (Left)", fsize=14, bold=True), column=0)
    
    tf_left = faces.TreeFace(trees[0], ts_left)
    tf_left.border.width = 0

    # --- 오른쪽 트리 설정 (180도 반전) ---
    ts_right = TreeStyle()
    ts_right.show_leaf_name = False
    ts_right.show_branch_support = False
    ts_right.scale = args.scale
    ts_right.branch_vertical_margin = 10
    ts_right.orientation = 1 # 💡 1 = 오른쪽에서 왼쪽으로 자람 (마주보기 핵심!)
    ts_right.title.add_face(TextFace("Tree 2 (Right) ", fsize=14, bold=True), column=0)
    
    tf_right = faces.TreeFace(trees[1], ts_right)
    tf_right.border.width = 0

    # 5. 가상의 마스터(Master) 트리를 만들고 양옆에 두 그림을 부착
    master = Tree(";") # 점 하나짜리 투명 트리
    ns_master = NodeStyle()
    ns_master["size"] = 0
    master.set_style(ns_master)

    # 왼쪽 그림 부착
    master.add_face(tf_left, column=0, position="branch-right")
    # 중앙에 빈 공간(Gap) 만들기
    master.add_face(TextFace("          "), column=1, position="branch-right") 
    # 오른쪽 그림 부착
    master.add_face(tf_right, column=2, position="branch-right")

    # 전체 캔버스 설정
    ts_master = TreeStyle()
    ts_master.show_leaf_name = False
    ts_master.show_scale = False
    ts_master.title.add_face(TextFace(" Face-to-Face Tree Comparison", fsize=18, bold=True), column=0)

    # 렌더링
    master.render(args.output, tree_style=ts_master, w=1400, units="px", dpi=150)
    print(f"🎉 두 개의 트리가 완벽하게 마주보는 시각화 완료: {args.output}")

if __name__ == "__main__":
    main()