import copy
from itertools import combinations, permutations, product
from collections import defaultdict
from tree import Node, isomorphic, generate_trees

# Generate all partitions in maximum <k> subsets
def generate_partitions(data, k):
    if len(data) == 1:
        yield [ data ]
        return
    first = data[0]
    for smaller in generate_partitions(data[1:], k):
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset] + smaller[n+1:]
        if len(smaller) < k:
            yield [ [ first ] ] + smaller

# Return all partitions of <n> elements into minimum <k_min> and maximum <k_max> subsets, sorted by number of subset
def get_sorted_parts(n, k_min, k_max):
    d = defaultdict(list)
    for p in generate_partitions([f"m{i}" for i in range(n)], k_max):
        if len(p) >= k_min:
            d[len(p)].append(p)
    return d

# Generate all possible tumor trees from <n_min> to <n_max> nodes with <k> mutations
def generate_tumors(n_min, n_max, k):
    if n_min > k:
        print(f"WARNING: trees starting at n_min={n_min} nodes can not contains only k={k} mutations")
    elif n_max > k:
        print(f"WARNING: trees ending at n_max={n_max} nodes can not contains only k={k} mutations")

    mutations_parts = get_sorted_parts(k, n_min, n_max)

    # For each possible number of node <n>
    for n, tree_batch in zip(range(n_min, n_max + 1), generate_trees(n_min, n_max)):
        # For each possible tree of <n> nodes for each partition of mutation in <n> subset
        for tree, label_part  in product(tree_batch, mutations_parts[n]):
            # For each possible node -> partition bijection
            for labels in permutations(label_part):
                # Generate the corresponding tree by assigning mutations labels
                new_tree = tree.deep_copy()
                for node, l in zip(new_tree.descendants, labels):
                    node.set_content(l)
                yield new_tree
            
          
K = 3
N = (1, 3)
print(f"starting generation with {K} mutations, from {N[0]} to {N[1]} nodes")
for t in generate_tumors(N[0], N[1], K):
    print("\nNew tree :")
    t.print_content()
    t.print_tree()
