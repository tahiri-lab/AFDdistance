import json
import re
from itertools import product
from typing import Counter, cast

import mp3treesim as mp3
import numpy as np
from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, get_all_mutations, precompute_all
from afd.utils import newick_to_tree
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from treelib.node import Node

from gmd.utils.distance_measures.AD_dist.ad import ancestor_descendant
from stereodist.CASet import caset, caset_intersection, caset_union
from stereodist.DISC import disc, disc_intersection, disc_union
from stereodist.utils import ancestor_sets
from utils import save_fig_safe

NORM_MUTS = None  # None, gene or position
NORM_FS = False  # for norm_muts potition only
PRINT_DEF = False
PRINT_NORM = False
PRINT_TREES = False
NORM_ROOT = True


with open("scripts/aml_trees.json", "r") as f:
    tree_infos = json.load(f)


def parse_mutation(mutation_id):
    # Regular expression for parsing mutation
    pattern = r"([-\+]?)([A-Za-z0-9]+)_p\.([A-Za-z])(\d+)([A-Za-z]+|fs)"
    match = re.match(pattern, mutation_id)

    if match:
        # Parse values
        loss = match.group(1)
        gene = match.group(2)
        initial_aa = match.group(3)
        position = int(match.group(4))
        new_aa = match.group(5)

        return {
            "loss": loss,
            "gene": gene,
            "initial_aa": initial_aa,
            "position": position,
            "new_aa": new_aa,
        }
    else:
        return None


def position_based(mutation_id, normalize_fs) -> str:
    data = parse_mutation(mutation_id)
    if data is None:
        return mutation_id
    if data["new_aa"] == "fs" and normalize_fs:
        data["position"] = "fs"
    return data["loss"] + data["gene"] + "_" + str(data["position"])


def gene_based(mutation_id) -> str:
    data = parse_mutation(mutation_id)
    if data is None:
        return mutation_id
    return data["loss"] + data["gene"]


def show_muts_content(trees_list) -> None:
    all_muts = dict()
    for t in trees_list:
        for n in t.all_nodes_itr():
            for m in n.data.muts:
                if not m in all_muts:
                    all_muts[m] = 0
                all_muts[m] += 1
    for m in sorted(all_muts.keys()):
        print(f"{all_muts[m]}\t{m}")


for t in tree_infos:
    # print(t["id"])
    new_tree = newick_to_tree(t["nwk"], freq=True, mode="bracket")
    if NORM_ROOT:
        r = cast(Node, new_tree[cast(str, new_tree.root)])
        r.data.freq = 0
    t["tree"] = new_tree

if PRINT_DEF:
    print("DEfault Trees :")
    show_muts_content([t["tree"] for t in tree_infos])


if NORM_MUTS is not None:
    for t in tree_infos:
        for n in t["tree"].all_nodes_itr():
            for i in range(len(n.data.muts)):
                if NORM_MUTS == "position":
                    n.data.muts[i] = position_based(n.data.muts[i], NORM_FS)
                elif NORM_MUTS == "gene":
                    n.data.muts[i] = gene_based(n.data.muts[i])
                else:
                    raise Exception("wrong norm_muts")

    if PRINT_NORM:
        print("\nNORM Trees :")
        show_muts_content([t["tree"] for t in tree_infos])


for t in tree_infos:
    t["tumortree"] = TumorTree(t["tree"])

if PRINT_TREES:
    for t in tree_infos:
        print(t["id"])
        t["tree"].show(data_property="muts")


all_trees = [t["tumortree"] for t in tree_infos]
# precompute_all(all_trees, [m for m in get_all_mutations(all_trees) if m != "Root"])
precompute_all(all_trees, get_all_mutations(all_trees))
maxdist_glob = afd_upper_bound(all_trees)


objects = np.array(all_trees)


def distanceafd(obj1, obj2):
    if afd_upper_bound([obj1, obj2], "val_rel_add") > 0:
        return afd_distance(obj1, obj2, False) / afd_upper_bound(
            [obj1, obj2], "val_rel_add"
        )
    else:
        return 1


def distancead(obj1, obj2):
    return ancestor_descendant(
        obj1.to_newick(freq=False, mode="bracket"),
        obj2.to_newick(freq=False, mode="bracket"),
    )


def distancecas(obj1, obj2):
    t1_ancestors = ancestor_sets(obj1.to_newick(freq=False, mode="bracket"))
    t2_ancestors = ancestor_sets(obj2.to_newick(freq=False, mode="bracket"))
    union = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))
    # union = [i for i in union if i != "Root"]

    for u in union:
        if u not in t1_ancestors:
            t1_ancestors[u] = set()
        if u not in t2_ancestors:
            t2_ancestors[u] = set()

    return caset(union, t1_ancestors, t2_ancestors)
    # return caset_union(
    #    obj1.to_newick(freq=False, mode="bracket"),
    #    obj2.to_newick(freq=False, mode="bracket"),
    # )


def distancecasi(obj1, obj2):
    t1_ancestors = ancestor_sets(obj1.to_newick(freq=False, mode="bracket"))
    t2_ancestors = ancestor_sets(obj2.to_newick(freq=False, mode="bracket"))
    intersection = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))
    # intersection = [i for i in intersection if i != "Root"]
    try:
        return caset(intersection, t1_ancestors, t2_ancestors)
        # return caset_intersection(
        #    obj1.to_newick(freq=False, mode="bracket"),
        #    obj2.to_newick(freq=False, mode="bracket"),
        # )
    except:
        return 1


def distancedisc(obj1, obj2):
    t1_ancestors = ancestor_sets(obj1.to_newick(freq=False, mode="bracket"))
    t2_ancestors = ancestor_sets(obj2.to_newick(freq=False, mode="bracket"))
    union = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))
    # union = [i for i in union if i != "Root"]

    for u in union:
        if u not in t1_ancestors:
            t1_ancestors[u] = set()
        if u not in t2_ancestors:
            t2_ancestors[u] = set()

    return disc(union, t1_ancestors, t2_ancestors)
    # return disc_union(
    #    obj1.to_newick(freq=False, mode="bracket"),
    #    obj2.to_newick(freq=False, mode="bracket"),
    # )


def distancedisci(obj1, obj2):
    t1_ancestors = ancestor_sets(obj1.to_newick(freq=False, mode="bracket"))
    t2_ancestors = ancestor_sets(obj2.to_newick(freq=False, mode="bracket"))
    intersection = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))
    # intersection = [i for i in intersection if i != "Root"]
    try:
        return disc(intersection, t1_ancestors, t2_ancestors)
        # return disc_intersection(
        #    obj1.to_newick(freq=False, mode="bracket"),
        #    obj2.to_newick(freq=False, mode="bracket"),
        # )
    except:
        return 1


def distancemp3(obj1, obj2):
    t1 = mp3.read_dotstring(obj1.to_gv())
    t2 = mp3.read_dotstring(obj2.to_gv())
    try:
        return 1 - mp3.similarity(t1, t2)
    except Exception:
        n1 = set(m for n in obj1.tree.all_nodes_itr() for m in n.data.muts)
        n2 = set(m for n in obj2.tree.all_nodes_itr() for m in n.data.muts)
        if len(n1.intersection(n2)) > 2:
            raise Exception
        return 1


metrics = {
    # "AFD": distanceafd,
    # "AD": distancead,
    # "CASet_U": distancecas,
    # "DISC_U": distancedisc,
    # "CASet_I": distancecasi,
    # "DISC_I": distancedisci,
    "MP3": distancemp3,
}


n = len(objects)
labels = np.array([t["id"] for t in tree_infos])

results: dict = {"labels": labels.tolist()}

for name, metric in metrics.items():

    results_metric = {}

    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = metric(objects[i], objects[j])
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    results_metric["dist"] = (
        dist_matrix.tolist()
    )  # [[float(x) for x in row] for row in dist_matrix.tolist()]

    condensed_dist = squareform(dist_matrix)
    linked = linkage(condensed_dist, method="average")

    # Test every possible cluster number

    results_metric["clusters"] = []

    for k in range(2, n):
        cluster_labels = fcluster(linked, k, criterion="maxclust")
        if len(set(cluster_labels)) == k:
            score = silhouette_score(dist_matrix, cluster_labels, metric="precomputed")
            inter = []
            intra = []
            for c1 in range(1, k + 1):
                indices_c1 = [i for i, x in enumerate(cluster_labels) if x == c1]
                for c2 in range(c1, k + 1):
                    indices_c2 = [i for i, x in enumerate(cluster_labels) if x == c2]
                    if c1 == c2:
                        intra += [
                            dist_matrix[i, j]
                            for i, j in product(indices_c1, indices_c2)
                            if i != j
                        ]
                    else:
                        inter += [
                            dist_matrix[i, j]
                            for i, j in product(indices_c1, indices_c2)
                            if i != j
                        ]

            results_metric["clusters"].append(
                {
                    "k": k,
                    "tags": cluster_labels.tolist(),
                    "sscore": score,
                    "inter": inter,
                    "intra": intra,
                }
            )

    results[name] = results_metric

    max_s = max(c["sscore"] for c in results_metric["clusters"])
    tags = next(c["tags"] for c in results_metric["clusters"] if c["sscore"] == max_s)

    tag_counts = Counter(tags)

    for t in tag_counts.keys():
        if tag_counts[t] <= 1:
            continue
        print("======================================")
        print("======================================")
        tagtrees = [objects[i] for i in range(len(objects)) if tags[i] == t]
        taglabels = [labels[i] for i in range(len(labels)) if tags[i] == t]
        for tt, tl in zip(tagtrees, taglabels):
            print(f"{tl} : ")
            tt.tree.show(data_property="displaydata")


# with open("scripts/aml-results-all.json", "w") as f:
# json.dump(results, f, indent=2)
