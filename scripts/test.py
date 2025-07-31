from itertools import product

from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree

from stereodist.CASet import caset_intersection, caset_union


def show_fd(tree: TumorTree, mutmap: dict[str, int]):
    for m1, m2 in product(mutmap.keys(), mutmap.keys()):
        if tree.fd[mutmap[m1], mutmap[m2]] != 0:
            print(f"FD({m1}, {m2}) = {tree.fd[mutmap[m1], mutmap[m2]]}")


t1 = "A:1(B:2(C:3(D:2,E:6),F:5(G:1)))"
t2 = "A:1(B:2(C:3(D:2,E:6),F:3(G:3)))"
t3 = "A:1(B:2(C:3(D:2,E:6),F:5(H:1)))"

tl1, tl2, tl3 = (newick_to_tree(t, freq=True, mode="bracket") for t in [t1, t2, t3])

tree1 = TumorTree.from_newick(t1, freq=True, mode="bracket")
tree2 = TumorTree.from_newick(t2, freq=True, mode="bracket")
tree3 = TumorTree.from_newick(t3, freq=True, mode="bracket")
mut_map = precompute_all([tree1, tree2, tree3])

print("===== T1 ======")
tl1.show(data_property="muts")
tl1.show(data_property="freq")
print("===== T2 ======")
tl2.show(data_property="muts")
tl2.show(data_property="freq")
print("===== T3 ======")
tl3.show(data_property="muts")
tl3.show(data_property="freq")

# Expect 0.538596491228
print(
    f"AFD(t1, t2) = {afd_distance(tree1, tree2)}/{afd_upper_bound([tree1, tree2])} = {afd_distance(tree1, tree2)/afd_upper_bound([tree1, tree2])}"
)

# Expect 6
print(
    f"AFD(t1, t3) = {afd_distance(tree1, tree3)}/{afd_upper_bound([tree1, tree3])} = {afd_distance(tree1, tree3)/afd_upper_bound([tree1, tree3])}"
)

# Expect 6
print(
    f"AFD(t2, t3) = {afd_distance(tree2, tree3)}/{afd_upper_bound([tree2, tree3])} = {afd_distance(tree2, tree3)/afd_upper_bound([tree2, tree3])}"
)


t4 = "A:1(B:1(X:1,C:1),D:1)"
t5 = "A:1(B:1(X:1,C:1),D:1(X:1))"  # duplicate x
t6 = "A:1(D:1(B:1(X:1,C:1)),D:1(X:1))"  # duplicate x and d

tl4, tl5, tl6 = (newick_to_tree(t, freq=True, mode="bracket") for t in [t4, t5, t6])

tree4 = TumorTree.from_newick(t4, freq=True, mode="bracket")
tree5 = TumorTree.from_newick(t5, freq=True, mode="bracket")
tree6 = TumorTree.from_newick(t6, freq=True, mode="bracket")
mut_map_2 = precompute_all([tree4, tree5, tree6])

print("===== T4 ======")
tl4.show(data_property="muts")
tl4.show(data_property="freq")
print("===== T5 ======")
tl5.show(data_property="muts")
tl5.show(data_property="freq")
print("===== T6 ======")
tl6.show(data_property="muts")
tl6.show(data_property="freq")


# expect 1.73
print(
    f"AFD(t4, t5) = {afd_distance(tree4, tree5)}/{afd_upper_bound([tree4, tree5])} = {afd_distance(tree4, tree5)/afd_upper_bound([tree4, tree5])}"
)

# expect 5.11
print(
    f"AFD(t4, t6) = {afd_distance(tree4, tree6)}/{afd_upper_bound([tree4, tree6])} = {afd_distance(tree4, tree6)/afd_upper_bound([tree4, tree6])}"
)

# expect 3.81
print(
    f"AFD(t5, t6) = {afd_distance(tree5, tree6)}/{afd_upper_bound([tree5, tree6])} = {afd_distance(tree5, tree6)/afd_upper_bound([tree5, tree6])}"
)
