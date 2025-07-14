from itertools import product

from scipy import sparse
from treelib.node import Node
from treelib.tree import Tree

from utils import SubCloneData, newick_to_tree


def compute_cumuled_freqs(tree: Tree) -> None:
    """Compute cumuled frequency on each node,
    The value is attached to node data in place.
    """
    for n_id in reversed(list(tree.expand_tree(mode=Tree.WIDTH))):
        node = tree[n_id]
        node.data.cumuled_freq = node.data.freq
        node.data.cumuled_freq += sum(
            [c.data.cumuled_freq for c in tree.children(n_id)], 0
        )


def node_frequency_diff(tree: Tree) -> sparse.lil_matrix:
    """Compute frequency difference rate between nodes,
    Value is returned as a sparse matrix.
    """
    compute_cumuled_freqs(tree)
    fd = sparse.lil_matrix((tree.size(), tree.size()))

    s = [str(tree.root)]
    marked = {n.identifier: False for n in tree.all_nodes_itr()}
    while len(s) > 0:
        # find the next child to compute
        n_id = s[-1]
        c_id = next(
            (c.identifier for c in tree.children(n_id) if not marked[c.identifier]),
            None,
        )

        # compute all frequency difference for c
        if c_id is not None:
            c_freq = tree[c_id].data.cumuled_freq
            # all ancestor of c or in the stack
            for a_id in s:
                a_freq = tree[a_id].data.cumuled_freq
                fd[int(a_id), int(c_id)] = c_freq / a_freq

            # mark c and put it on the stack to compute its children frequency diffs
            marked[c_id] = True
            s.append(c_id)

        # all children of n were already computed so discard it
        else:
            s.pop()

    return fd


def mut_frequency_diff(tree: Tree, muts_map: dict) -> sparse.lil_matrix:
    """Compute frequency difference rate between mutations
    Value is returned as a sparse matrix.
    muts_map map each mutation to an integer from 0 to m-1
    """
    nfd = node_frequency_diff(tree)
    m = len(muts_map)
    mfd = sparse.lil_matrix((m, m))

    for i, j in zip(*nfd.nonzero()):
        val = nfd[i, j]
        for m_i, m_j in product(tree[str(i)].data.muts, tree[str(j)].data.muts):
            idx_i, idx_j = muts_map[m_i], muts_map[m_j]
            mfd[idx_i, idx_j] = val

    return mfd


def fd_diff(fd1: sparse.lil_matrix, fd2: sparse.lil_matrix) -> float:
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


def afd(input1: Tree | str, input2: Tree | str):
    if type(input1) == str:
        tree1 = newick_to_tree(input1)
    else:
        tree1 = Tree(input1)

    if type(input2) == str:
        tree2 = newick_to_tree(input2)
    else:
        tree2 = Tree(input2)

    all_muts = sum([n.data.muts for n in tree1.all_nodes_itr()], [])
    all_muts += sum([n.data.muts for n in tree2.all_nodes_itr()], [])
    muts_map = {m: i for i, m in enumerate(set(all_muts))}
    m = len(muts_map)

    fd1 = mut_frequency_diff(tree1, muts_map)
    fd2 = mut_frequency_diff(tree2, muts_map)

    return fd_diff(fd1, fd2) / (m * (m - 1) / 2)
