#!/usr/bin/env python3
import os
import sys
import argparse

# Render on a headless server without GUI if needed
os.environ['QT_QPA_PLATFORM'] = 'offscreen' 

from ete3 import Tree, TreeStyle, NodeStyle

def main():
    # 1. 인자 파서(Argument Parser) 설정
    parser = argparse.ArgumentParser(description="Root and render a phylogenetic tree using ETE3.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input Newick tree file (.nwk or .tre).")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for the output SVG file path (e.g., /path/to/output/GeneA).")
    parser.add_argument("-g", "--outgroup", required=True, help="String to identify the outgroup tips (e.g., 'Sphagnum' or 'Bryum_argenteum').")
    
    args = parser.parse_args()

    # 2. 트리 파일 로드
    try:
        with open(args.input, "r") as f:
            newick_str = f.read().strip()
        t = Tree(newick_str)
        print(f"✔ Successfully loaded tree from: {args.input}")
    except Exception as e:
        print(f"❌ Error loading tree: {e}")
        sys.exit(1)

    # 3. Outgroup 팁 식별 (단일 종 또는 클레이드)
    outgroup_nodes = []
    for leaf in t.get_leaves():
        if args.outgroup in leaf.name:
            outgroup_nodes.append(leaf)

    if not outgroup_nodes:
        print(f"❌ Error: No tips containing '{args.outgroup}' were found in the tree.")
        sys.exit(1)

    # 4. 트리 루팅(Rooting)
    try:
        if len(outgroup_nodes) == 1:
            # 단일 종만 매칭될 경우 해당 노드로 루팅
            t.set_outgroup(outgroup_nodes[0])
        else:
            # 여러 종이 매칭될 경우 이들의 가장 가까운 공통 조상(LCA)을 찾아 루팅
            lca = t.get_common_ancestor(outgroup_nodes)
            t.set_outgroup(lca)
        
        print(f"✔ Successfully rooted tree using {len(outgroup_nodes)} node(s) matching '{args.outgroup}'.")
    except Exception as e:
        print(f"❌ Error during rooting: {e}")
        sys.exit(1)

    # 5. 시각화 스타일 설정
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = True 
    ts.branch_vertical_margin = 10 
    
    # 커스텀 노드 스타일
    nstyle = NodeStyle()
    nstyle["size"] = 0 
    nstyle["vt_line_width"] = 2 
    nstyle["hz_line_width"] = 2 
    nstyle["vt_line_color"] = "#333333" 
    nstyle["hz_line_color"] = "#333333"

    for n in t.traverse():
        n.set_style(nstyle)

    # 6. 렌더링 및 파일 저장
    output_file = f"{args.output_prefix}_rooted_tree.svg"
    
    # 출력할 디렉토리가 없다면 생성
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        t.render(output_file, w=183, units="mm", tree_style=ts)
        print(f"✔ Tree rendered and saved to: {output_file}")
    except Exception as e:
        print(f"❌ Error rendering tree: {e}")

if __name__ == "__main__":
    main()