import random
import sys
from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree, tree_to_newick

import distances as dist
from stereodist.CASet import caset_union
from stereodist.DISC import disc_union
from utils import save_fig_safe

random.seed(sys.argv[1])

# parse base tree
nwk = "(((((((((J)I)H)G)F)E)D)C)B)A"
# nwk = "(((((((((J:1)I:2)H:3)G:4)F:5)E:6)D:7)C:8)B:9)A:10"
tbase = newick_to_tree(nwk, freq=False, mode="bracket", cst="rnd")
# tbase = newick_to_tree(nwk, freq=True, mode="bracket", cst="1")

# create reference
ref_tree = TumorTree(deepcopy(tbase))
ref_nwk = tree_to_newick(tbase, freq=False, mode="bracket")

# generate tree variations
nwk_variants = []
tree_variants = []
for i in range(0, 9):
    newtree = deepcopy(tbase)
    newtree.move_node("9", str(i))
    nwk_variants.append(tree_to_newick(newtree, freq=False, mode="bracket"))
    tree_variants.append(TumorTree(newtree))

# precompute fd matrices
precompute_all([ref_tree] + tree_variants)

# base dataframe
data = pd.DataFrame(
    data={
        "var_id": list(range(len(tree_variants))),
    }
)

tree_variants = list(reversed(tree_variants))

# metrics
data["AFD"] = [dist.distanceafd(ref_tree, v) for v in tree_variants]
data["CASet"] = [dist.distancecas(ref_tree, v) for v in tree_variants]
data["DISC"] = [dist.distancedisc(ref_tree, v) for v in tree_variants]
data["AD"] = [dist.distancead(ref_tree, v) for v in tree_variants]
data["MP3"] = [dist.distancemp3(ref_tree, v) for v in tree_variants]

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
plt.xticks(list(range(0, len(tree_variants))))
plt.xlabel("Position")

# plt.show()
save_fig_safe("figures/depthtest.pdf")

data.to_csv("lesfilsdechineng.csv")
# tbase = deepcopy(base_trunk)
# tbase.create_node()
