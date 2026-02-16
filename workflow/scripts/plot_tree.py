import sys
from ete3 import Tree, TreeStyle, TextFace, faces

input_tree = snakemake.input[0]
output_img = snakemake.output[0]
sample_name = snakemake.wildcards.sample

def get_custom_ts():
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.title.add_face(TextFace(f"Phylogeny: {sample_name}", fsize=12), column=0)
    return ts

def layout(node):
    if node.is_leaf():
        name_face = TextFace(node.name, fsize=10)
        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    else:
        if node.support:
            color = "red" if node.support < 80 else "blue"
            support_face = TextFace(f"{node.support:.0f}", fsize=8, fgcolor=color)
            faces.add_face_to_node(support_face, node, column=0, position="branch-top")

try:
    t = Tree(input_tree, format=1)

    try:
        t.set_outgroup(t.get_midpoint_outgroup())
    except:
        pass # pass if error occurs or if tree is too small

    ts = get_custom_ts()
    ts.layout_fn = layout

    # save image
    t.render(output_img, tree_style=ts, w=800, units="px")

except Exception as e:
    print(f"Error plotting {sample_name}: {e}")
    sys.exit(1)