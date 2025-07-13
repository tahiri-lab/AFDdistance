from afd import afd
from utils import newick_to_tree

t1 = "A:1(B:2(C:3(D:2,E:6),F:5(G:1)))"
t2 = "A:1(B:2(C:3(D:2,E:6),F:3(G:3)))"

tree1 = newick_to_tree(t1)
tree2 = newick_to_tree(t2)

tree1.show(data_property="muts")
print("T1:")
tree1.show(data_property="freq")
print("T2:")
tree2.show(data_property="freq")

# Expect 0.538596491228
print(f"AFD(t1, t2) = {afd(t1, t2)}")
