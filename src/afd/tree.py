"""
Useful things to play around with tumor trees
"""

from collections import defaultdict
from itertools import product
from typing import Counter, Iterable

from scipy import sparse
from treelib.node import Node
from treelib.tree import Tree

from .utils import SubCloneData, newick_to_tree, tree_to_gv, tree_to_newick


class TumorTree:
    _tree: Tree
    _mutations_map: dict[str, list[int]]

    _cumuled_freqs: dict[str, float]
    _fd_local: sparse.lil_matrix
    _fd_global: sparse.lil_matrix | None = None
    _global_map: dict[str, int]

    def __init__(self, tree: Tree):
        if isinstance(tree, str):
            self._tree = newick_to_tree(tree)
        else:
            self._tree = tree

        self.remap_mutations()
        self.compute_cumuled_freq()
        self.compute_local_fd()
        self.unmap_mutations()

    @property
    def mutations(self) -> list[str]:
        """List of mutations of the tree"""
        return list(self._mutations_map.keys())

    @property
    def tree(self) -> Tree:
        """Tree object of the tree"""
        return self._tree

    @property
    def fd(self) -> sparse.lil_matrix[float]:
        """Return the global fd matrix"""
        if self._fd_global == None:
            raise Exception(
                "Accessing uncomputed FD matrix, compute fd_global prior to accessing it"
            )
        return self._fd_global

    @property
    def size(self) -> int:
        """size in total number of mutations (including duplicates)"""
        return len(sum((ids for ids in self._mutations_map.values()), []))

    def count_mutation(self, mut: str) -> int:
        """return the number of occurence of a mutation in the tree"""
        return len(self._mutations_map.get(mut, []))

    def count_relation(self, mut1: str, mut2: str) -> int:
        """return the number of occurences of a relation in the tree"""
        return len(
            [
                1
                for m1, m2 in product(
                    self._mutations_map[mut1], self._mutations_map[mut2]
                )
                if self._fd_local[m1, m2] > 0
            ]
        )

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

    def unmap_mutations(self) -> None:
        for n in self.tree.all_nodes_itr():
            for i, m in enumerate(n.data.muts):
                m_name = next(
                    (k for k, v in self._mutations_map.items() if m in v), "none"
                )
                n.data.muts[i] = m_name

    def compute_cumuled_freq(self) -> dict:
        self._cumuled_freqs = dict()
        for nid in reversed(list(self.tree.expand_tree(mode=Tree.WIDTH))):
            node = self.tree[nid]
            self._cumuled_freqs[nid] = node.data.freq
            for child in self.tree.children(nid):
                cid = child.identifier
                self._cumuled_freqs[nid] += self._cumuled_freqs[cid]
            node.data.cumuled_freq = self._cumuled_freqs[nid]
        return self._cumuled_freqs

    def compute_local_fd(self) -> sparse.lil_matrix:
        """Compute the FD relation matrix local to the tree.
        which mean mutation indetifier are mapped only locally, and
        each mutation occurence has a unique id
        This is roughly in $O(m^2)$ because we fill one cell of the matrix
        for each ancestor descendant relationship. The set of ancestor descendant
        relationship upper bound is $O(m^2)$
        """
        m = self.size
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
                for mut2 in node_muts:
                    if mut2 != mut:
                        self._fd_local[mut, mut2] = 1

        return self._fd_local

    def merge_recurrent(self, m1, m2) -> float:
        """Compute the final falue for ancestry m1 to m2 by
        collapsing every value that relate to this relation.
        This gives a sort of sub-matrix with each row relating to m1 and
        each column relating to m2.
        - First each row is collapsed (relation from the same ancestor (row))
        - Then the collapsed values are summed (relations from different ancestors m1)
        """
        total = 0.0

        for m1_row in self._mutations_map[m1]:
            all_values = [
                self._fd_local[m1_row, m2_col]
                for m2_col in self._mutations_map[m2]
                if self._fd_local[m1_row, m2_col] != 0
            ]
            collapsed_value = len(all_values) * sum(all_values, 0)
            total += collapsed_value

        return total

    def compute_global_fd(self, global_map: dict[str, int]) -> sparse.lil_matrix:
        """This map the local fd matrix to a global matrix which represent relationship
        for any mutation being in the global set of trees to compare.
        """
        self._global_map = global_map
        m = len(global_map)
        self._fd_global = sparse.lil_matrix((m, m))

        for mut_i, mut_j in product(self.mutations, self.mutations):
            if mut_i not in global_map.keys() or mut_j not in global_map.keys():
                continue
            i_glob = global_map[mut_i]
            j_glob = global_map[mut_j]

            self._fd_global[i_glob, j_glob] = self.merge_recurrent(mut_i, mut_j)

        return self._fd_global

    def to_newick(self, freq: bool, mode: str):
        """Compute the newick string for the tree
        Args:
            freq: add frequency for each subclone as <subclone>**:0.1**
            mode: wether to group multi-labels in "brackets" or to separate
                them by space ("spaced")
        """
        return tree_to_newick(self.tree, freq, mode)

    def to_gv(self) -> str:
        return tree_to_gv(self.tree)

    @classmethod
    def from_newick(cls, nwk: str, freq: bool, mode: str, cst: str = "1"):
        return cls(newick_to_tree(nwk, freq, mode, cst))


def get_all_mutations(trees: Iterable[TumorTree]) -> list[str]:
    """Return the union of all mutations in a set of trees"""
    muts = set()
    for t in trees:
        muts.update(t.mutations)
    return list(muts)


def precompute_all(
    trees: Iterable[TumorTree], mutations: list | None = None
) -> dict[str, int]:
    """For every tree, compute global map with respect to the total set of mutations
    so they can be compared together flawlessly
    """
    if mutations is None:
        mutations = get_all_mutations(trees)

    global_map = {m: i for i, m in enumerate(mutations)}
    for t in trees:
        t.compute_global_fd(global_map)
    return global_map
