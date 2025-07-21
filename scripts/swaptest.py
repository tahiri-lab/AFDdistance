from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from treelib.tree import Tree

from afd import afd
from stereodist.CASet import caset_intersection, caset_union
from stereodist.DISC import disc_intersection, disc_union
from utils import newick_to_tree, tree_to_newick


def swap_labels(tree: Tree, nid1: str, nid2: str) -> Tree:
    newtree = deepcopy(tree)
    n1 = newtree.get_node(nid1)
    n2 = newtree.get_node(nid2)
    if n1 is None or n2 is None:
        raise Exception()
    temp = n1.data.muts[0]
    n1.data.muts[0] = n2.data.muts[0]
    n2.data.muts[0] = temp
    return newtree


tbase = newick_to_tree("((C:1,D:1)B:1,(F:1,G:1)E:1)A:1")
newick_tbase = tree_to_newick(tbase)

swap1 = [
    swap_labels(tbase, "5", "6"),
    swap_labels(tbase, "2", "3"),
]
swap2 = [
    swap_labels(tbase, "1", "2"),
    swap_labels(tbase, "1", "3"),
    swap_labels(tbase, "4", "5"),
    swap_labels(tbase, "4", "6"),
]
swap3 = [
    swap_labels(tbase, "2", "5"),
    swap_labels(tbase, "2", "6"),
    swap_labels(tbase, "3", "5"),
    swap_labels(tbase, "3", "6"),
]
swap4 = [
    swap_labels(tbase, "1", "5"),
    swap_labels(tbase, "1", "6"),
    swap_labels(tbase, "4", "2"),
    swap_labels(tbase, "4", "3"),
]
swap5 = [
    swap_labels(tbase, "1", "4"),
]
swap6 = [
    swap_labels(tbase, "0", "1"),
    swap_labels(tbase, "0", "4"),
]
swap7 = [
    swap_labels(tbase, "0", "2"),
    swap_labels(tbase, "0", "3"),
    swap_labels(tbase, "0", "5"),
    swap_labels(tbase, "0", "6"),
]

swaps = [swap1, swap2, swap3, swap4, swap5, swap6, swap7]
swaps_flat = sum(swaps, [])
newick_swaps_flat = [tree_to_newick(s) for s in swaps_flat]

swap_classes = sum([[f"swap{i}"] * len(s) for i, s in enumerate(swaps)], [])
swaps_ids = list(range(len(swaps_flat)))

data = pd.DataFrame(
    data={
        "swap_id": swaps_ids,
        "swap_class": swap_classes,
    }
)

data["afd"] = [afd(tbase, s) for s in swaps_flat]
data["caset_u"] = [caset_union(newick_tbase, s) for s in newick_swaps_flat]
data["disc_u"] = [disc_union(newick_tbase, s) for s in newick_swaps_flat]

df_melted = pd.melt(
    data,
    id_vars=["swap_id", "swap_class"],
    value_vars=["afd", "caset_u", "disc_u"],
    var_name="distance metric",
    value_name="distance value",
)

fig, ax = plt.subplots(figsize=(10, 5))

sns.lineplot(
    data=df_melted,
    x="swap_id",
    y="distance value",
    hue="distance metric",
    linestyle="--",
    marker="o",
)
# plt.show()
plt.savefig("../figures/swap.pdf")
