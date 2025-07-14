from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from treelib.tree import Tree

from afd import afd
from stereodist.CASet import caset_union
from stereodist.DISC import disc_union
from utils import newick_to_tree, tree_to_newick

tbase = newick_to_tree("(((((((((J:1)I:1)H:1)G:1)F:1)E:1)D:1)C:1)B:1)A:1")
newick_tbase = tree_to_newick(tbase)

variants = []
for i in range(0, 9):
    newtree = deepcopy(tbase)
    newtree.move_node("9", str(i))
    variants.append(newtree)
newick_variants = [tree_to_newick(v) for v in variants]

var_ids = list(range(len(variants)))

data = pd.DataFrame(
    data={
        "var_id": var_ids,
    }
)

data["afd"] = [afd(tbase, v) for v in variants]
data["caset_u"] = [caset_union(newick_tbase, v) for v in newick_variants]
data["disc_u"] = [disc_union(newick_tbase, v) for v in newick_variants]

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
plt.show()


# tbase = deepcopy(base_trunk)
# tbase.create_node()
