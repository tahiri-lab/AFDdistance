from itertools import product

from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree

with open("scripts/data/fimo_crc2_newick.txt", "r") as f:
    fimo = f.readline()
with open("scripts/data/scarlet_crc2_newick.txt", "r") as f:
    scarlet = f.readline()
with open("scripts/data/sciclonefit_crc2_newick.txt", "r") as f:
    scf = f.readline()

fimo_tree = TumorTree.from_newick(fimo, False, "space")
scarlet_tree = TumorTree.from_newick(scarlet, False, "space")
scf_tree = TumorTree.from_newick(scf, False, "space")

trees = [fimo_tree, scarlet_tree, scf_tree]
trees_name = ["fimo", "scarlet", "sciclonefit"]

precompute_all(trees)
maxdist_glob = afd_upper_bound(trees)

for i in range(len(trees)):
    for j in range(i + 1, len(trees)):
        print(f"AFD({trees_name[i]}, {trees_name[j]}) = ", end="")
        print(
            f"{round(afd_distance(trees[i], trees[j]), 3)}/{round(afd_upper_bound([trees[i], trees[j]]), 3)} = {round(afd_distance(trees[i], trees[j])/afd_upper_bound([trees[i], trees[j]]), 3)}"
        )
        print(f"AFD({trees_name[i]}, {trees_name[j]}) = ", end="")
        print(
            f"{round(afd_distance(trees[i], trees[j]), 3)}/{maxdist_glob} = {round(afd_distance(trees[i], trees[j])/maxdist_glob, 3)}"
        )
