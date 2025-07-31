from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree, read_newicks, tree_to_newick

from stereodist.CASet import caset_intersection, caset_union
from stereodist.DISC import disc_intersection, disc_union

gain_1_nwk = read_newicks("scripts/data/8mut1gainsA.txt")
gain_2_nwk = read_newicks("scripts/data/8mut2gainsA.txt")
gain_3_nwk = read_newicks("scripts/data/8mut3gainsA.txt")
gain_4_nwk = read_newicks("scripts/data/8mut4gainsA.txt")

gain_1_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket") for n in gain_1_nwk
]
gain_2_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket") for n in gain_2_nwk
]
gain_3_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket") for n in gain_3_nwk
]
gain_4_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket") for n in gain_4_nwk
]
all_trees = gain_1_trees + gain_2_trees + gain_3_trees + gain_4_trees
precompute_all(all_trees)

# Compute with local normalization
data = []
for t1, t2 in product(gain_1_trees, gain_1_trees):
    data.append({"ngain": 1, "afd": afd_distance(t1, t2) / afd_upper_bound([t1, t2])})
for t1, t2 in product(gain_1_trees, gain_2_trees):
    data.append({"ngain": 2, "afd": afd_distance(t1, t2) / afd_upper_bound([t1, t2])})
for t1, t2 in product(gain_1_trees, gain_3_trees):
    data.append({"ngain": 3, "afd": afd_distance(t1, t2) / afd_upper_bound([t1, t2])})
for t1, t2 in product(gain_1_trees, gain_4_trees):
    data.append({"ngain": 4, "afd": afd_distance(t1, t2) / afd_upper_bound([t1, t2])})

df = pd.DataFrame(data=data)
sns.lineplot(data=df, x="ngain", y="afd")
plt.title("AFD with normalization for each pair")
plt.show()


# Compute with global normalization
global_upper_bound = afd_upper_bound(all_trees)
data2 = []
for t1, t2 in product(gain_1_trees, gain_1_trees):
    data2.append({"ngain": 1, "afd": afd_distance(t1, t2) / global_upper_bound})
for t1, t2 in product(gain_1_trees, gain_2_trees):
    data2.append({"ngain": 2, "afd": afd_distance(t1, t2) / global_upper_bound})
for t1, t2 in product(gain_1_trees, gain_3_trees):
    data2.append({"ngain": 3, "afd": afd_distance(t1, t2) / global_upper_bound})
for t1, t2 in product(gain_1_trees, gain_4_trees):
    data2.append({"ngain": 4, "afd": afd_distance(t1, t2) / global_upper_bound})


df2 = pd.DataFrame(data=data2)
sns.lineplot(data=df2, x="ngain", y="afd")
plt.title("AFD with set normalization")
plt.show()


all_dists = np.zeros((len(all_trees), len(all_trees)))
for i in range(len(all_trees)):
    for j in range(len(all_trees)):
        dist = afd_distance(all_trees[i], all_trees[j])
        all_dists[i, j] = dist / global_upper_bound
        all_dists[j, i] = dist / global_upper_bound

sns.heatmap(all_dists, cmap="coolwarm", annot=False, cbar=True)
plt.show()
