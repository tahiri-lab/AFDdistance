import json
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from gmd.utils.distance_measures.AD_dist.ad import ancestor_descendant
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from stereodist.CASet import caset_intersection, caset_union
from stereodist.DISC import disc_intersection, disc_union

from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree, precompute_all
from afd.utils import newick_to_tree

NORM_MUTS = None  # None, gene or position
NORM_FS = True  # for norm_muts potition only
PRINT_DEF = False
PRINT_NORM = True
PRINT_TREES = False


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
precompute_all(all_trees)
maxdist_glob = afd_upper_bound(all_trees)


objects = np.array(all_trees)


def distanceafd(obj1, obj2):
    return afd_distance(obj1, obj2, True) / afd_upper_bound([obj1, obj2], "nrel")


def distanceafd2(obj1, obj2):
    return afd_distance(obj1, obj2, False) / afd_upper_bound([obj1, obj2], "sumrel")


def distancead(obj1, obj2):
    return ancestor_descendant(
        obj1.to_newick(freq=False, mode="bracket"),
        obj2.to_newick(freq=False, mode="bracket"),
    )


def distancecas(obj1, obj2):
    return caset_union(
        obj1.to_newick(freq=False, mode="bracket"),
        obj2.to_newick(freq=False, mode="bracket"),
    )


def distancedisc(obj1, obj2):
    return disc_union(
        obj1.to_newick(freq=False, mode="bracket"),
        obj2.to_newick(freq=False, mode="bracket"),
    )


metrics = {
    "AFD max 1 (nrel)": distanceafd,
    "AFD no max (sumrel)": distanceafd2,
    "AD": distancead,
    "CASet U": distancecas,
    "DISC U": distancedisc,
}

for name, metric in metrics.items():
    print(f"\n================= {name} =================\n")
    n = len(objects)
    labels = np.array([t["id"] for t in tree_infos])

    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = metric(objects[i], objects[j])
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    condensed_dist = squareform(dist_matrix)
    linked = linkage(condensed_dist, method="average")

    # === Silhouette Score vs Number of Clusters ===
    max_clusters = n - 1  # Avoid clusters > n
    scores = []

    for k in range(2, max_clusters + 1):
        cluster_labels = fcluster(linked, k, criterion="maxclust")
        score = silhouette_score(dist_matrix, cluster_labels, metric="precomputed")
        scores.append(score)

    n_clusters = scores.index(max(scores)) + 2

    # Get cluster assignments
    clusters = fcluster(linked, n_clusters, criterion="maxclust")

    # Create DataFrame with distance matrix
    df = pd.DataFrame(dist_matrix, index=labels, columns=labels)

    # Sort by cluster assignments
    cluster_order = pd.DataFrame({"label": labels, "cluster": clusters}).sort_values(
        "cluster"
    )
    sorted_labels = cluster_order["label"].tolist()
    df_sorted = df.loc[sorted_labels, sorted_labels]
    # Create row colors based on clusters

    colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
    cluster_colors = [
        colors[cluster - 1] for cluster in cluster_order["cluster"]
    ]  # Create clustermap

    sns.clustermap(
        df_sorted,
        row_cluster=False,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        row_colors=cluster_colors,
        col_colors=cluster_colors,
    )
    plt.title(f"{name}: CLustering map")

    for cluster_id in sorted(set(clusters)):
        cluster_labels = [labels[i] for i, c in enumerate(clusters) if c == cluster_id]
        print(f"Cluster {cluster_id}: {cluster_labels}")

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
    dendrogram(
        linked,
        labels=labels,  # Optional: replace with more meaningful labels
        distance_sort=True,
        show_leaf_counts=True,
        ax=ax1,
    )
    ax1.set_title(f"{name}: Hierarchical Clustering Dendrogram")
    ax1.set_xlabel("Object Index")
    ax1.set_ylabel("Distance")
    # plt.show()

    # Plot silhouette score
    ax2.plot(range(2, max_clusters + 1), scores, marker="o")
    ax2.set_title(f"{name}: Silhouette Score vs Number of Clusters")
    ax2.set_xlabel("Number of Clusters")
    ax2.set_ylabel("Silhouette Score")
    ax2.grid(True)
    ax2.set_xticks(range(2, max_clusters + 1))

    plt.show()

#
# from fcmeans import FCM
# from scipy.spatial.distance import squareform
# from sklearn.manifold import MDS
#
# metrics = {
#     "AFD max 1 (nrel)": distanceafd,
#     "AFD no max (sumrel)": distanceafd2,
# }
#
# for name, metric in metrics.items():
#     n = len(objects)
#     labels = [t["id"] for t in tree_infos]
#
#     dist_matrix = np.zeros((n, n))
#     for i in range(n):
#         for j in range(i + 1, n):
#             d = metric(objects[i], objects[j])
#             dist_matrix[i, j] = d
#             dist_matrix[j, i] = d
#
#     condensed_dist = squareform(dist_matrix)
#     linked = linkage(condensed_dist, method="average")
#
#     min_clusters = 2
#     max_clusters = n - 1
#     scores = []
#     clusterings = []
#
#     for k in range(min_clusters, max_clusters + 1):
#         fcm = FCM(
#             m=2.0, max_iter=150, error=1e-5, n_jobs=1, n_clusters=k, random_state=0
#         )
#         fcm.fit(dist_matrix)  # It expects features, so we’ll transform below
#
#         # Use hard labels from fuzzy membership for silhouette
#         labels_hard = fcm.u.argmax(axis=1)
#
#         # Silhouette requires a distance matrix, use precomputed mode
#         score = silhouette_score(dist_matrix, labels_hard, metric="precomputed")
#         scores.append(score)
#         clusterings.append(fcm)
#
#     best_index = np.argmax(scores)
#     best_k = best_index + min_clusters
#     best_model = clusterings[best_index]
#     best_labels = best_model.u.argmax(axis=1)
#
#     fig, axes = plt.subplots(
#         1, 3, figsize=(18, 6), gridspec_kw={"width_ratios": [1, 1, 1.2]}
#     )
#     ax1, ax2, ax3 = axes
#
#     # Plot 1: Silhouette Score
#     ax1.plot(range(min_clusters, max_clusters + 1), scores, marker="o")
#     ax1.set_title("Silhouette Score vs Number of Clusters")
#     ax1.set_xlabel("Number of Clusters")
#     ax1.set_ylabel("Silhouette Score")
#     ax1.grid(True)
#     ax1.set_xticks(range(min_clusters, max_clusters + 1))
#
#     # Plot 2: MDS projection (2D embedding)
#     mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0)
#     coords = mds.fit_transform(dist_matrix)
#     scatter = ax2.scatter(coords[:, 0], coords[:, 1], c=best_labels, cmap="tab10", s=50)
#     ax2.set_title(f"MDS Projection (k={best_k})")
#     ax2.set_xticks([])
#     ax2.set_yticks([])
#     legend_labels = np.unique(best_labels)
#     ax2.legend(*scatter.legend_elements(), title="Cluster")
#
#     # Plot 3: Heatmap of reordered distance matrix
#     # Order objects by cluster assignment
#     sort_idx = np.argsort(best_labels)
#     sorted_matrix = dist_matrix[sort_idx][:, sort_idx]
#     sorted_labels = [labels[i] for i in sort_idx] if labels else sort_idx
#
#     sns.heatmap(
#         sorted_matrix,
#         ax=ax3,
#         cmap="viridis",
#         xticklabels=sorted_labels,
#         yticklabels=sorted_labels,
#     )
#     ax3.set_title("Distance Matrix (Cluster-Sorted)")
#     ax3.tick_params(labelrotation=90)
#
#     plt.tight_layout()
#     plt.show()
