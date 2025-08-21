from .tree import TumorTree


def afd_upper_bound(trees: list[TumorTree]) -> float:
    """Compute the maximum distance between two trees
    Maximum distance is the sum of all FD values in both trees
    This function assume already precomputed trees
    """
    return sum((float(t.fd.sum()) for t in trees), 0.0)


def afd_distance(
    t1: TumorTree, t2: TumorTree
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
        diff += abs(val1 - val2)

    return diff
