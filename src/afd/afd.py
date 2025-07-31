from scipy import sparse

from .tree import TumorTree, get_all_mutations


def afd_upper_bound(trees: list[TumorTree], type="nmuts") -> float:
    """Given the number of mutation m there are m*(m-1)/2 ancestor descendant relations
    with maximum difference of 1.
    For each mutations we count its number of occurence in the tree where it is the most frequent
    """
    if type == "nmuts":
        all_muts = get_all_mutations(trees)
        m = len(all_muts)
        return m * (m - 1)
    elif type == "nrel":
        return sum((t.fd.nnz for t in trees), 0)
    elif type == "sumrel":
        return sum((t.fd.sum() for t in trees), 0)
    else:
        return 0


def afd_distance(
    t1: TumorTree, t2: TumorTree, valmax: bool = True
):  # , mutations: dict[str, int]) -> float:
    """Assume t1 and t2 with already precomputed fd_global matrices"""
    fd1, fd2 = t1.fd, t2.fd

    diff = 0

    coords1 = set(zip(*fd1.nonzero()))
    coords2 = set(zip(*fd2.nonzero()))
    all_coords = coords1.union(coords2)

    for i, j in sorted(all_coords):
        val1 = float(fd1[i, j])
        val2 = float(fd2[i, j])

        # relation not in t1, count the number of occurences in t2
        if valmax and val1 == 0:
            m1 = next((k for k, v in t2._global_map.items() if v == i), "")
            m2 = next((k for k, v in t2._global_map.items() if v == j), "")
            diff += t2.count_relation(m1, m2)
            # diff += 1
        # relation not in t2, count the number of occurences in t1
        elif valmax and val2 == 0:
            m1 = next((k for k, v in t1._global_map.items() if v == i), "")
            m2 = next((k for k, v in t1._global_map.items() if v == j), "")
            diff += t1.count_relation(m1, m2)
            # diff += 1
        else:
            diff += abs(val1 - val2)

    return diff
