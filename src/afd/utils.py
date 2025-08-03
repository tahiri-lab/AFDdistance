"""Handle translation between newick format and treelib trees
In newick the content (mutations) of a given node are either:
- "space": A B C D
- "bracket": {A,B,C,D}
Trees can be read / write from both format with or without frequency data.
Frequency data is represented by a float such as: {A,B,C}:0.83
"""

import re
from typing import cast

from treelib.node import Node
from treelib.tree import Tree


class SubCloneData:
    """Struct to store the node data in tumor tree"""

    def __init__(self, freq: float, mutation_load: list) -> None:
        self.freq = freq
        self.cumuled_freq = None
        self.muts = mutation_load

    @property
    def displaydata(self) -> str:
        return f"{self.muts} = {self.cumuled_freq}"


def parent_children_split(newick: str) -> tuple[str, str | None]:
    """Split parent(children) parts of the newick string
    If children part is None, node is a leaf
    """
    if newick.startswith("("):
        last_index = newick.rfind(")")
        return newick[last_index + 1 :], newick[1:last_index]
    if newick.endswith(")"):
        first_index = newick.find("(")
        return newick[:first_index:], newick[first_index + 1 : -1]
    else:
        return newick, None


def get_mutation_load(
    subclone: str, parsefreq: bool, mode: str
) -> tuple[float, list[str]]:
    """Get the list of mutation in a subclone from newick"""
    if parsefreq:
        muts, freq = subclone.split(":")
    else:
        muts, freq = subclone, 1

    if mode == "bracket":
        match = re.search(r"\{(.*?)\}", muts)
        if match:
            return float(freq), match.group(1).split(",")
        else:
            return float(freq), [muts]
    elif mode == "space":
        mutlist = [m.strip() for m in muts.split(" ") if m.strip() != ""]
        return float(freq), mutlist
    else:
        raise Exception(f"Unknown subclone mode {mode}")


# https://bitbucket.org/oesperlab/stereodist/src/master/utils.py
def outer_comma_split(newick):
    """Split every child.
    Split a Newick subtring on commas, ignoring those contained in parentheses and brackets.
    :param newick: the string to split
    """
    chunk_start = 0
    parens = 0
    brackets = 0

    for i, char in enumerate(newick):
        if char == "(":
            parens += 1
        elif char == ")":
            parens -= 1
        elif char == "{":
            brackets += 1
        elif char == "}":
            brackets -= 1
        elif char == "," and parens == 0 and brackets == 0:
            yield newick[chunk_start:i]
            chunk_start = i + 1

    yield newick[chunk_start:]


def parse_newick(
    newick: str,
    tree: Tree,
    parent: Node | None,
    parsefreq: bool,
    mode: str = "bracket",
) -> None:
    """Recursively parse the newick string to fill the tree.
    fill <tree> in place
    """
    node, childrens = parent_children_split(newick)
    freq, muts = get_mutation_load(node, parsefreq, mode)
    p_node = tree.create_node(
        str(tree.size()),
        str(tree.size()),
        parent=parent,
        data=SubCloneData(freq, muts),
    )

    if childrens is not None:
        for child in outer_comma_split(childrens):
            parse_newick(child, tree, p_node, parsefreq, mode)


def newick_to_tree(newick: str, freq: bool = True, mode: str = "bracket") -> Tree:
    """Parse a newick string to create a treelib Tree
    return the tree with root "0"
    """
    t = Tree()
    parse_newick(newick, t, None, freq, mode)
    t.root = "0"
    return t


def tree_to_newick(
    tree: Tree, freq: bool = False, mode: str = "bracket", node: Node | None = None
) -> str:
    """Recursively parse a tree to create its newick string"""
    if node is None:
        node = cast(Node, tree.get_node(tree.root))

    if len(node.data.muts) > 1:
        if mode == "bracket":
            parent = "{" + ",".join(map(str, node.data.muts)) + "}"
        elif mode == "space":
            parent = " ".join(map(str, node.data.muts))
        else:
            raise Exception()
    else:
        parent = str(node.data.muts[0])

    if freq:
        parent += f":{node.data.freq}"

    childrens = []
    for c in tree.children(node.identifier):
        childrens.append(tree_to_newick(tree, freq, mode, c))

    result = ""
    if len(childrens) > 0:
        result += "(" + ",".join(childrens) + ")"
    result += parent

    return result


def read_newicks(file: str) -> list[str]:
    trees = []
    with open(file, mode="r") as f:
        for line in f:
            nwk = line.strip().strip(";")
            if nwk != "":
                trees.append(nwk)
    return trees


def write_newicks(newicks: list[str], file: str):
    with open(file, mode="w") as f:
        f.writelines([f"{nwk};\n" for nwk in newicks])


def convert(file, source, dest):
    """convert a newick file from space to bracket or vis versa"""
    trees = [newick_to_tree(t, False, source) for t in read_newicks(file)]
    nwks = [tree_to_newick(t, False, dest) for t in trees]
    new_file = ".".join(file.split(".")[:-1]) + f"_{dest}." + file.split(".")[-1]
    write_newicks(nwks, new_file)
