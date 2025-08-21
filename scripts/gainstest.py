"""Replicate of the test from the GMD article"""

# TODO: add results for gmd adjusted methods

from itertools import product

import distances as dist
import numpy as np
import pandas as pd
import seaborn as sns
from utils import save_fig_safe

from afd.tree import TumorTree, get_all_mutations, precompute_all
from afd.utils import read_newicks

gain_1_nwk = read_newicks("scripts/data/8mut1gainsA.txt")
gain_2_nwk = read_newicks("scripts/data/8mut2gainsA.txt")
gain_3_nwk = read_newicks("scripts/data/8mut3gainsA.txt")
gain_4_nwk = read_newicks("scripts/data/8mut4gainsA.txt")

gain_1_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket", cst="rnd") for n in gain_1_nwk
]
gain_2_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket", cst="rnd") for n in gain_2_nwk
]
gain_3_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket", cst="rnd") for n in gain_3_nwk
]
gain_4_trees = [
    TumorTree.from_newick(n, freq=False, mode="bracket", cst="rnd") for n in gain_4_nwk
]
all_trees = gain_1_trees + gain_2_trees + gain_3_trees + gain_4_trees
precompute_all(all_trees, get_all_mutations(all_trees))

data = []
for t1, t2 in product(gain_1_trees, gain_1_trees):
    data.append(
        {
            "ngain": 1,
            "AFD": dist.distanceafd(t1, t2),
            "CASet": dist.distancecas(t1, t2),
            "DISC": dist.distancedisc(t1, t2),
            "AD": dist.distancead(t1, t2),
            "MP3": dist.distancemp3(t1, t2),
        }
    )
for t1, t2 in product(gain_1_trees, gain_2_trees):
    data.append(
        {
            "ngain": 2,
            "AFD": dist.distanceafd(t1, t2),
            "CASet": dist.distancecas(t1, t2),
            "DISC": dist.distancedisc(t1, t2),
            "AD": dist.distancead(t1, t2),
            "MP3": dist.distancemp3(t1, t2),
        }
    )
for t1, t2 in product(gain_1_trees, gain_3_trees):
    data.append(
        {
            "ngain": 3,
            "AFD": dist.distanceafd(t1, t2),
            "CASet": dist.distancecas(t1, t2),
            "DISC": dist.distancedisc(t1, t2),
            "AD": dist.distancead(t1, t2),
            "MP3": dist.distancemp3(t1, t2),
        }
    )
for t1, t2 in product(gain_1_trees, gain_4_trees):
    data.append(
        {
            "ngain": 4,
            "AFD": dist.distanceafd(t1, t2),
            "CASet": dist.distancecas(t1, t2),
            "DISC": dist.distancedisc(t1, t2),
            "AD": dist.distancead(t1, t2),
            "MP3": dist.distancemp3(t1, t2),
        }
    )


df = pd.DataFrame(data=data)

df_melted = pd.melt(
    df,
    id_vars=["ngain"],
    value_vars=["AFD", "CASet", "DISC", "AD", "MP3"],
    var_name="distance method",
    value_name="distance value",
)

sns.lineplot(data=df_melted, x="ngain", y="distance value", hue="distance method")
save_fig_safe("results/simulated/gaintest_plot.pdf")


all_dists = np.zeros((len(all_trees), len(all_trees)))
for i in range(len(all_trees)):
    for j in range(len(all_trees)):
        distij = dist.distanceafd(all_trees[i], all_trees[j])
        all_dists[i, j] = distij
        all_dists[j, i] = distij

sns.heatmap(all_dists, cmap="coolwarm", annot=False, cbar=True)
save_fig_safe("results/simulated/gaintest_hm_afd.pdf")
