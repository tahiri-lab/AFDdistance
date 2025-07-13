import re
from typing import cast

from treelib.node import Node
from treelib.tree import Tree


class SubCloneData:
    def __init__(self, freq: int, mutation_load: list) -> None:
        self.freq = freq
        self.cumuled_freq = None
        self.muts = mutation_load


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


def get_mutation_load(subclone: str) -> list[str]:
    """Get the list of mutation in a subclone from newick"""
    brackets = r"\{(.*?)\}"
    match = re.search(brackets, subclone)
    if match:
        return match.group(1).split(",")
    else:
        return [subclone]


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


def parse_newick(newick: str, tree: Tree, parent: Node | None) -> None:
    """Recursively parse the newick string to fill the tree.
    fill <tree> in place
    """
    node, childrens = parent_children_split(newick)
    p_node = tree.create_node(
        str(tree.size()),
        str(tree.size()),
        parent=parent,
        data=SubCloneData(0, get_mutation_load(node)),
    )

    if childrens is not None:
        for child in outer_comma_split(childrens):
            parse_newick(child, tree, p_node)


def newick_to_tree(newick: str) -> Tree:
    """Parse a newick string to create a treelib Tree
    return the tree with root "0"
    """
    t = Tree()
    parse_newick(newick, t, None)
    t.root = "0"
    return t


def tree_to_newick(tree: Tree, node: Node | None = None) -> str:
    """Recursively parse a tree to create its newick string"""
    if node is None:
        node = cast(Node, tree.get_node(tree.root))

    if len(node.data.muts) > 1:
        parent = "{" + ",".join(map(str, node.data.muts)) + "}"
    else:
        parent = str(node.data.muts[0])

    childrens = []
    for c in tree.children(node.identifier):
        childrens.append(tree_to_newick(tree, c))

    result = ""
    if len(childrens) > 0:
        result += "(" + ",".join(childrens) + ")"
    result += parent

    return result


# t = newick_to_tree("(MAP2K1,(LRRC16A){PLA2G16,EXOC6B}){GPR158,SLC12A1}")
# t.show(data_property="muts")
# t.show()
#
# print(tree_to_newick(t))
