from scipy import sparse

from .tree import TumorTree, get_all_mutations

# def fd_matrix(tree: Tree, map: dict[str, list[int]]) -> sparse.lil_matrix:
#     m = len(sum(map.values(), []))
#     fd = sparse.lil_matrix((m, m))
#
#     return fd
#
#
# def node_frequency_diff(tree: Tree) -> sparse.lil_matrix:
#     """Compute frequency difference rate between nodes,
#     Value is returned as a sparse matrix.
#     """
#     compute_cumuled_freqs(tree)
#     fd = sparse.lil_matrix((tree.size(), tree.size()))
#
#     s = [str(tree.root)]
#     marked = {n.identifier: False for n in tree.all_nodes_itr()}
#     while len(s) > 0:
#         # find the next child to compute
#         n_id = s[-1]
#         c_id = next(
#             (c.identifier for c in tree.children(n_id) if not marked[c.identifier]),
#             None,
#         )
#
#         # compute all frequency difference for c
#         if c_id is not None:
#             c_freq = tree[c_id].data.cumuled_freq
#             # all ancestor of c or in the stack
#             for a_id in s:
#                 a_freq = tree[a_id].data.cumuled_freq
#                 fd[int(a_id), int(c_id)] += c_freq / a_freq
#
#             # mark c and put it on the stack to compute its children frequency diffs
#             marked[c_id] = True
#             s.append(c_id)
#
#         # all children of n were already computed so discard it
#         else:
#             s.pop()
#
#     return fd
#
#
# def mut_frequency_diff(tree: Tree, muts_map: dict) -> sparse.lil_matrix:
#     """Compute frequency difference rate between mutations
#     Value is returned as a sparse matrix.
#     muts_map map each mutation to an integer from 0 to m-1
#     """
#     nfd = node_frequency_diff(tree)
#     m = len(muts_map)
#     mfd = sparse.lil_matrix((m, m))
#
#     for i, j in zip(*nfd.nonzero()):
#         val = nfd[i, j]
#         for m_i, m_j in product(tree[str(i)].data.muts, tree[str(j)].data.muts):
#             idx_i, idx_j = muts_map[m_i], muts_map[m_j]
#             mfd[idx_i, idx_j] = val
#
#     return mfd
#


# TODO: handle maximum value in recurrent mutations context
def fd_diff(fd1: sparse.lil_matrix, fd2: sparse.lil_matrix) -> float:
    """FD diff account for the difference between each fd value in both matrix
    if a fd value is zero in one matrix but not the other, maximum value is assigned.
    Without recurrent mutation the maximum value is 1 but this is to be reconsidered.
    """
    diff = 0

    coords1 = set(zip(*fd1.nonzero()))
    coords2 = set(zip(*fd2.nonzero()))
    all_coords = coords1.union(coords2)

    for i, j in sorted(all_coords):
        val1 = float(fd1[i, j])
        val2 = float(fd2[i, j])
        if val1 == 0 or val2 == 0:
            diff += 1.0
        else:
            diff += abs(val1 - val2)

    return diff


def afd_upper_bound(t1: TumorTree, t2: TumorTree) -> float:
    """Given the number of mutation m there are m*(m-1)/2 ancestor descendant relations
    with maximum difference of 1.
    Each duplicate is counted as a new mutation
    """
    m = len(get_all_mutations([t1, t2])) + t1.n_duplicate + t2.n_duplicate
    return m * (m - 1) / 2


def afd_distance(t1: TumorTree, t2: TumorTree) -> float:
    """Assume t1 and t2 with already precomputed fd_global matrices"""

    return fd_diff(t1.fd, t2.fd)
