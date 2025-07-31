from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from stereodist.CASet import caset_union
from stereodist.DISC import disc_union
from utils import save_fig_safe

from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree, tree_to_newick

# parse base tree
nwk = "(((((((((J:1)I:1)H:1)G:1)F:1)E:1)D:1)C:1)B:1)A:1"
tbase = newick_to_tree(nwk, freq=True, mode="bracket")

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

# metrics
data["afd"] = [
    afd_distance(ref_tree, v) / afd_upper_bound([ref_tree, v]) for v in tree_variants
]
data["caset_u"] = [caset_union(ref_nwk, v) for v in nwk_variants]
data["disc_u"] = [disc_union(ref_nwk, v) for v in nwk_variants]

df_melted = pd.melt(
    data,
    id_vars=["var_id"],
    value_vars=["afd", "caset_u", "disc_u"],
    var_name="distance metric",
    value_name="distance value",
)

sns.lineplot(
    data=df_melted,
    x="var_id",
    y="distance value",
    hue="distance metric",
    linestyle="--",
    marker="o",
)

# plt.show()
save_fig_safe("figures/depthtest.pdf")


# tbase = deepcopy(base_trunk)
# tbase.create_node()
