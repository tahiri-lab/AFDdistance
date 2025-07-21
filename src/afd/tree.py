"""
Useful things to play around with tumor trees
"""

from collections import defaultdict
from itertools import product
from typing import Counter, Iterable

from scipy import sparse
from treelib.node import Node
from treelib.tree import Tree

from .utils import SubCloneData, newick_to_tree, tree_to_newick


class TumorTree:
    _tree: Tree
    _mutations_map: dict[str, list[int]]

    _cumuled_freqs: dict[str, float]
    _fd_local: sparse.lil_matrix
    _fd_global: sparse.lil_matrix | None = None

    def __init__(self, tree: Tree):
        if isinstance(tree, str):
            self._tree = newick_to_tree(tree)
        else:
            self._tree = tree

        self.remap_mutations()
        self.compute_cumuled_freq()
        self.compute_local_fd()

    @property
    def mutations(self) -> list[str]:
        """List of mutations of the tree"""
        return list(self._mutations_map.keys())

    @property
    def tree(self) -> Tree:
        """Tree object of the tree"""
        return self._tree

    @property
    def fd(self) -> sparse.lil_matrix:
        if self._fd_global == None:
            raise Exception(
                "Accessing uncomputed FD matrix, compute fd_global prior to accessing it"
            )
        return self._fd_global

    @property
    def n_duplicate(self) -> int:
        """count the number of duplicate for each mutation
        without counting the first occurence"""
        return sum((len(i) - 1 for i in self._mutations_map.values()), 0)

    def remap_mutations(self) -> None:
        """Map each mutation to an integer from 0 to m
        and rename them on the tree.

        If a mutation appear multiple times, each occurence
        is mapped onto a new integer.
        """
        identifier = 0
        self._mutations_map = defaultdict(list)
        for n in self.tree.all_nodes_itr():
            for i, m in enumerate(n.data.muts):
                self._mutations_map[m].append(identifier)
                n.data.muts[i] = identifier
                identifier += 1

    def compute_cumuled_freq(self) -> dict:
        self._cumuled_freqs = dict()
        for nid in reversed(list(self.tree.expand_tree(mode=Tree.WIDTH))):
            node = self.tree[nid]
            self._cumuled_freqs[nid] = node.data.freq
            for child in self.tree.children(nid):
                cid = child.identifier
                self._cumuled_freqs[nid] += self._cumuled_freqs[cid]
        return self._cumuled_freqs

    def compute_local_fd(self) -> sparse.lil_matrix:
        """Compute the FD relation matrix local to the tree.
        which mean mutation indetifier are mapped only locally, and
        each mutation occurence has a unique id
        This is roughly in $O(m^2)$ because we fill one cell of the matrix
        for each ancestor descendant relationship. The set of ancestor descendant
        relationship upper bound is $O(m^2)$
        """
        m = self._tree.size()
        self._fd_local = sparse.lil_matrix((m, m))

        # walk in reverse so we can infer parent descendants from children descendants
        for nid in reversed(list(self.tree.expand_tree(mode=Tree.WIDTH))):
            node = self.tree[nid]
            node_muts = node.data.muts

            # FD row for the mutations in the current node
            # (they all have the same set of descendants)
            node_desc_rates = sparse.lil_matrix((1, m))

            for child in self.tree.children(nid):
                cid = child.identifier
                child_muts = child.data.muts

                # Get node to child rate and add FD value for direct descendant mutations
                node_child_rate = self._cumuled_freqs[cid] / self._cumuled_freqs[nid]
                for mut in child_muts:
                    node_desc_rates[0, mut] = node_child_rate

                # Add FD value for indirect descendant mutations
                node_desc_rates += node_child_rate * self._fd_local[child_muts[0]]

            # write the FD values for each mutation of the current node
            for mut in node_muts:
                self._fd_local[mut] = node_desc_rates

        # collapse columns (same relation from the same ancestor)
        # final result is assigned to the first id of each unique mutation
        for cols in self._mutations_map.values():
            if len(cols) == 1:
                continue
            for row in self._fd_local.rows:
                all_values = [
                    self._fd_local[row, c] for c in cols if self._fd_local[row, c] != 0
                ]
                collapsed_value = len(all_values) * sum(all_values, 0)
                self._fd_local[row, cols[0]] = collapsed_value
                for c in cols[1:]:
                    self._fd_local[row, c] = 0

        # collaps rows (same relations from different ancestors)
        # final result is assigned to the first id of each unique mutation
        for rows in self._mutations_map.values():
            if len(rows) == 1:
                continue
            for r in rows[1:]:
                self._fd_local[rows[0]] += self._fd_local[r]
                self._fd_local[r] = 0

        return self._fd_local

    def compute_global_fd(self, global_map: dict[str, int]) -> sparse.lil_matrix:
        """This map the local fd matrix to a global matrix which represent relationship
        for any mutation being in the global set of trees to compare.
        """
        m = len(global_map)
        self._fd_global = sparse.lil_matrix((m, m))

        for mut_i, mut_j in product(self.mutations, self.mutations):
            if mut_i not in global_map.keys() or mut_j not in global_map.keys():
                continue
            i_loc = self._mutations_map[mut_i][0]
            j_loc = self._mutations_map[mut_j][0]
            i_glob = global_map[mut_i]
            j_glob = global_map[mut_j]

            self._fd_global[i_glob, j_glob] = self._fd_local[i_loc, j_loc]

        return self._fd_global

    def to_newick(self, freq: bool, mode: str):
        """Compute the newick string for the tree
        Args:
            freq: add frequency for each subclone as <subclone>**:0.1**
            mode: wether to group multi-labels in "brackets" or to separate
                them by space ("spaced")
        """
        return tree_to_newick(self.tree, freq, mode)

    @classmethod
    def from_newick(cls, nwk: str, freq: bool, mode: str):
        return cls(newick_to_tree(nwk, freq, mode))


def get_all_mutations(trees: Iterable[TumorTree]) -> list[str]:
    """Return the union of all mutations in a set of trees"""
    muts = set()
    for t in trees:
        muts.update(t.mutations)
    return list(muts)


def precompute_all(trees: Iterable[TumorTree]) -> dict[str, int]:
    """For every tree, compute global map with respect to the total set of mutations
    so they can be compared together flawlessly
    """
    global_map = {m: i for i, m in enumerate(get_all_mutations(trees))}
    for t in trees:
        t.compute_global_fd(global_map)
    return global_map
