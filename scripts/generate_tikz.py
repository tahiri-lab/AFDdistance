import json
import os
import sys
from copy import deepcopy

import numpy as np
from afd.utils import newick_to_tree
from treelib.node import Node
from treelib.tree import Tree

colors_base = {
    "red!50": 0,
    "blue!50": 0,
    "green!50": 0,
    "yellow!50": 0,
    "orange!50": 0,
    "pink!50": 0,
    "brown!50": 0,
    "gray!50": 0,
    "black!50": 0,
    "cyan!50": 0,
    "purple!50": 0,
    "magenta!50": 0,
    "lime!50": 0,
    "violet!50": 0,
}

colors = deepcopy(colors_base)


def get_new_color():
    mink, minv = "red!50", colors["red!50"]
    for k, v in colors.items():
        if v < minv:
            mink, minv = k, v
    colors[mink] += 1
    return mink


def get_tikz(tree: Tree, node: Node, level: int, norm: int) -> str:
    if node.data.muts[0] == "Root":
        res = f"node [treenode] (root) {{root}}"
    else:
        res = f"node [treenode] {{{node.data.muts[0].replace("_", "\\_")}}}"
    for child in tree.children(node.identifier):
        res += (
            "\n"
            + (level + 1) * "\t"
            + "child {"
            + get_tikz(tree, child, level + 1, norm)
            + "}"
        )
    size = np.log2(node.data.freq / norm * 4.5 + 1)
    size_lab = int(node.data.freq)
    color = get_new_color()
    if node.data.muts[0] != "Root":
        res += (
            "\n"
            + (level + 1) * "\t"
            + f"child {{node [sizenode={{{color}}}{{{size}cm}}] {{{size_lab}}}}}"
        )
        res += "\n" + "\t" * level
    return res


def gen_trees_tikz(tree_ids, cluster):
    global colors
    for patient in tree_ids:
        colors = deepcopy(colors_base)
        with open("scripts/aml_trees.json", mode="r") as f:
            data = json.load(f)

        newick = next((d["nwk"] for d in data if d["id"] == patient), None)
        if newick is None:
            raise Exception
        tree = newick_to_tree(newick, freq=True, mode="bracket")
        norm = max(
            n.data.freq for n in tree.all_nodes_itr() if n.data.muts[0] != "Root"
        )

        tikz = """
        \\documentclass[tikz]{standalone}
        \\usepackage{tikz}
        \\usetikzlibrary{trees,shapes.geometric,positioning}
        
        \\begin{document}
        \\begin{tikzpicture}[
            sizenode/.style 2 args={circle, draw=#1, fill=#1, minimum size=#2},
            treenode/.style={ellipse, draw, inner sep=2pt, align=center},
            level 1/.style={level distance=2cm,sibling distance=3cm},
            level 2/.style={level distance=2cm,sibling distance=3cm}
        ]
        
        """

        tikz += "\\" + get_tikz(tree, tree[tree.root], 0, norm) + ";\n\n"

        tikz += f"\\node [left=1cm of root] {{\\textbf{{{patient}}}}};\n\n"

        tikz += """
        \\end{tikzpicture}
        \\end{document}
        """

        os.makedirs(f"figures/trees/{cluster}", exist_ok=True)

        with open(f"figures/trees/{cluster}/tikz-{patient}.tex", "w") as f:
            f.write(tikz)


gen_trees_tikz(["AML-05-001", "AML-69-001"], "cluster0")
gen_trees_tikz(["AML-13-001", "AML-40-001", "AML-52-001", "AML-63-001"], "cluster1")
gen_trees_tikz(["AML-01-001", "AML-03-001", "AML-49-001"], "cluster2")
gen_trees_tikz(
    ["AML-09-001", "AML-12-001", "AML-14-001", "AML-31-001", "AML-74-001"], "cluster3"
)
gen_trees_tikz(["AML-06-001", "AML-35-001", "AML-53-001"], "cluster4")
gen_trees_tikz(["AML-24-001", "AML-65-001", "AML-68-001"], "cluster5")
gen_trees_tikz(["AML-19-001", "AML-43-001"], "cluster6")
gen_trees_tikz(["AML-17-001", "AML-37-001"], "cluster7")
gen_trees_tikz(["AML-20-001", "AML-27-001", "AML-28-001", "AML-51-001"], "cluster8")
gen_trees_tikz(["AML-29-001", "AML-60-001"], "cluster9")
gen_trees_tikz(["AML-26-001", "AML-30-001"], "cluster10")
gen_trees_tikz(["AML-08-001", "AML-58-001", "AML-64-001"], "cluster11")
gen_trees_tikz(["AML-07-001", "AML-33-001"], "cluster12")
