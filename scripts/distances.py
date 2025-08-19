import mp3treesim as mp3
from afd import afd_distance, afd_upper_bound

from gmd.utils.distance_measures.AD_dist.ad import ancestor_descendant
from stereodist.CASet import caset, caset_intersection, caset_union
from stereodist.DISC import disc, disc_intersection, disc_union


def distanceafd(obj1, obj2):
    return afd_distance(obj1, obj2, False) / afd_upper_bound(
        [obj1, obj2], "val_rel_add"
    )


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
