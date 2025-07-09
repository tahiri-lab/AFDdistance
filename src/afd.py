from copy import deepcopy

import numpy as np
from scipy import sparse

from tree import Node


def compute_fd(tree: Node, freqs: dict):
    """Return fd values as a sparse array
    - tree: the root of the tree, each node name must be an integer from
        0 to tree.size() - 1
    - freqs: frequency of each node in a dict, with key being node name.
    """
    freqs_cum = deepcopy(freqs)
    fd = sparse.lil_matrix((tree.size(), tree.size()))

    for n in tree.postorder():
        for c in n.children:
            freqs_cum[n.name] += freqs_cum[c.name]

    print("CUMULED FREQUENCIES\n")
    print(freqs_cum)
    print()

    s = [tree]
    marked = {n.name: False for n in tree.descendants}
    while len(s) > 0:
        n = s[-1]
        unmarked = next((c for c in n.children if not marked[c.name]), None)
        if unmarked:
            fd[n.name, unmarked.name] = freqs_cum[unmarked.name] / freqs_cum[n.name]
            for ancestor in s[:-1]:
                fd[ancestor.name, unmarked.name] = (
                    freqs_cum[unmarked.name] / freqs_cum[ancestor.name]
                )
            marked[unmarked.name] = True
            s.append(unmarked)
        else:
            s.pop()

    return fd


def afd(fd1, fd2):
    fd1_coo = fd1.tocoo()
    fd2_coo = fd2.tocoo()

    coords1 = set(zip(fd1_coo.row, fd1_coo.col))
    coords2 = set(zip(fd2_coo.row, fd2_coo.col))
    all_coords = coords1 | coords2

    afd = 0
    for row, col in all_coords:
        val1 = fd1[row, col]
        val2 = fd2[row, col]
        afd += abs(val1 - val2)


n0 = Node(0)
n1 = Node(1, parent=n0)
n2 = Node(2, parent=n1)
n3 = Node(3, parent=n2)
n4 = Node(4, parent=n2)
n5 = Node(5, parent=n1)
n6 = Node(6, parent=n5)

f = {0: 1, 1: 2, 2: 2, 3: 5, 4: 2, 5: 1, 6: 7}

print("INITIAL TREE\n")
n0.print_tree()
print()
print(n0.size())

print("PER NODE FREQUENCY\n")
print(f)
print()

fd = compute_fd(n0, f)
print("FREQUENCY DIFFERENCE\n")
for i in range(0, 7):
    for j in range(0, 7):
        if fd[i, j] != 0:
            print(f"FD({i},{j}) = {fd[i,j]}")
