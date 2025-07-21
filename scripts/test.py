from itertools import product

from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree


def show_fd(tree: TumorTree, mutmap: dict[str, int]):
    for m1, m2 in product(mutmap.keys(), mutmap.keys()):
        if tree.fd[mutmap[m1], mutmap[m2]] != 0:
            print(f"FD({m1}, {m2}) = {tree.fd[mutmap[m1], mutmap[m2]]}")


t1 = "A:1(B:2(C:3(D:2,E:6),F:5(G:1)))"
t2 = "A:1(B:2(C:3(D:2,E:6),F:3(G:3)))"
t3 = "A:1(B:2(C:3(D:2,E:6),F:5(H:1)))"

tree1 = TumorTree.from_newick(t1, freq=True, mode="bracket")
tree2 = TumorTree.from_newick(t2, freq=True, mode="bracket")
tree3 = TumorTree.from_newick(t3, freq=True, mode="bracket")
mut_map = precompute_all([tree1, tree2, tree3])

tree1.tree.show(data_property="muts")
print("T1:")
tree1.tree.show(data_property="freq")
print("T2:")
tree2.tree.show(data_property="freq")
print("===\n")
print("T3:")
tree3.tree.show(data_property="muts")
tree3.tree.show(data_property="freq")

# Expect 0.538596491228
print(
    f"AFD(t1, t2) = {afd_distance(tree1, tree2)}/{afd_upper_bound(tree1, tree2)} = {afd_distance(tree1, tree2)/afd_upper_bound(tree1, tree2)}"
)

# Expect 6
print(
    f"AFD(t1, t3) = {afd_distance(tree1, tree3)}/{afd_upper_bound(tree1, tree3)} = {afd_distance(tree1, tree3)/afd_upper_bound(tree1, tree3)}"
)

# Expect 6 too
print(
    f"AFD(t2, t3) = {afd_distance(tree2, tree3)}/{afd_upper_bound(tree2, tree3)} = {afd_distance(tree2, tree3)/afd_upper_bound(tree2, tree3)}"
)
