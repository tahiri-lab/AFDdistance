"""Replicate of the test from the CASet / DISC article
with AFD results.
should add a reference here probably but this will do it for now
"""

from copy import deepcopy

import distances as dist
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from utils import save_fig_safe

from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree

# parse base tree
nwk = "(((((((((J)I)H)G)F)E)D)C)B)A"
tbase = newick_to_tree(nwk, freq=False, mode="bracket", cst="1")

# generate tree variations
trees = []
for i in range(0, 9):
    newtree = deepcopy(tbase)
    newtree.move_node("9", str(i))
    trees.append(TumorTree(newtree))

trees = list(reversed(trees))
ref_tree = trees[0]

# precompute fd matrices
precompute_all(trees)

# base dataframe
data = pd.DataFrame(
    data={
        "var_id": list(range(len(trees))),
    }
)

# metrics
data["AFD"] = [dist.distanceafd(ref_tree, v) for v in trees]
data["CASet"] = [dist.distancecas(ref_tree, v) for v in trees]
data["DISC"] = [dist.distancedisc(ref_tree, v) for v in trees]
data["AD"] = [dist.distancead(ref_tree, v) for v in trees]
data["MP3"] = [dist.distancemp3(ref_tree, v) for v in trees]

df_melted = pd.melt(
    data,
    id_vars=["var_id"],
    value_vars=["AFD", "CASet", "DISC", "AD", "MP3"],
    var_name="Mesure",
    value_name="Distance",
)

sns.lineplot(
    data=df_melted,
    x="var_id",
    y="Distance",
    hue="Mesure",
    linestyle="--",
    marker="o",
)
plt.xticks(list(range(0, len(trees))))
plt.xlabel("Position")

# plt.show()
save_fig_safe("results/simulated/depthtest.pdf")
data.to_csv("results/simulated/depthtest.csv")
