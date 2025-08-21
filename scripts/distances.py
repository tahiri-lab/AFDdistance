"""Helpers so any distance call is straightforward and
all can be done from TumorTree object"""

import mp3treesim as mp3
from gmd.utils.distance_measures.AD_dist.ad import ancestor_descendant
from stereodist.CASet import caset, caset_intersection, caset_union
from stereodist.DISC import disc, disc_intersection, disc_union

from afd import afd_distance, afd_upper_bound
from afd.tree import TumorTree


def distanceafd(tree1: TumorTree, tree2: TumorTree):
    return afd_distance(tree1, tree2) / afd_upper_bound([tree1, tree2])


def distancead(tree1: TumorTree, tree2: TumorTree):
    return ancestor_descendant(
        tree1.to_newick(freq=False, mode="bracket"),
        tree2.to_newick(freq=False, mode="bracket"),
    )


def distancecas(tree1: TumorTree, tree2: TumorTree):
    return caset_union(
        tree1.to_newick(freq=False, mode="bracket"),
        tree2.to_newick(freq=False, mode="bracket"),
    )


def distancedisc(tree1: TumorTree, tree2: TumorTree):
    return disc_union(
        tree1.to_newick(freq=False, mode="bracket"),
        tree2.to_newick(freq=False, mode="bracket"),
    )


def distancemp3(tree1: TumorTree, tree2: TumorTree):
    t1 = mp3.read_dotstring(tree1.to_gv())
    t2 = mp3.read_dotstring(tree2.to_gv())
    try:
        return 1 - mp3.similarity(t1, t2)
    except Exception:
        n1 = set(m for n in tree1.tree.all_nodes_itr() for m in n.data.muts)
        n2 = set(m for n in tree2.tree.all_nodes_itr() for m in n.data.muts)
        if len(n1.intersection(n2)) > 2:
            raise Exception
        return 1
