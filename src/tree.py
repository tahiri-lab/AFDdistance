import copy
from itertools import combinations, permutations


class Node:
    def __init__(self, name, content=[], parent=None):
        self.name = name
        self.content = content
        self.parent = None
        self.children = []
        self.is_leaf = True
        self.descendants = [self]
        self.depth = parent.depth + 1 if parent else 0
        if parent:
            self.attach(parent)

    # update depth recursively
    def set_depth(self, d):
        self.depth = d
        for c in self.children:
            c.set_depth(d + 1)

    def get_node(self, name):
        return next((n for n in self.descendants if n.name == name), None)

    def size(self):
        return len(self.descendants)

    def add_descendants(self, desc):
        self.descendants += desc
        if self.parent is not None:
            self.parent.add_descendants(desc)

    def remove_descendants(self, desc):
        self.descendants = [d for d in self.descendants if not d in desc]
        if self.parent is not None:
            self.parent.remove_descendants(desc)

    def get_leaves(self):
        return [d for d in self.descendants if d.is_leaf]

    # prune the subtree rooted at self
    def prune(self):
        if self.parent is not None:
            self.parent.children.remove(self)
            self.parent.remove_descendants(self.descendants)
            self.parent.is_leaf = len(self.parent.children) > 0
            self.parent = None
        self.depth = 0

    # attach a parentless node to a new parent
    def attach(self, node):
        if self.parent is not None:
            raise Exception(
                "Can not attach a node that already has a parent. forgot to prune() ?"
            )
        self.parent = node
        self.parent.children.append(self)
        self.parent.add_descendants(self.descendants)
        self.parent.is_leaf = False
        self.set_depth(node.depth + 1)

    # prune the subtree rooted at <self> and regrapht it to node
    def spr(self, node):
        self.prune()
        self.attach(node)

    def add_content(self, content):
        self.content += content

    def set_content(self, content):
        self.content = content

    def get_neighbors(self):
        return self.children + [self.parent] if self.parent else []

    def deep_copy(self):
        return copy.deepcopy(self)

    def print_content(self):
        print("Contents :")
        for d in self.descendants:
            print(f"\t{d.name} : {d.content}")

    def print_tree(self, level=0):
        indend = "  " * level
        print(f"{indend}{self.name}")
        for c in self.children:
            c.print_tree(level + 1)

    def postorder(self):
        f = [self]
        s = []
        while len(f) > 0:
            current_node = f.pop()
            s.append(current_node)
            for c in current_node.children:
                f.append(c)
        s.reverse()
        return s


# Tell wether two trees are isomorphic
def isomorphic(T1, T2):
    if len(T1.children) != len(T1.children):
        return False
    if len(T1.descendants) != len(T2.descendants):
        return False

    n_child = len(T1.children)
    if n_child == 0:
        return True

    for t1_children in permutations(T1.children):
        all_match = True
        for i in range(n_child):
            if not isomorphic(t1_children[i], T2.children[i]):
                all_match = False
                break
        if all_match:
            return True

    return False


# Generate all trees shape from k_min to k_max nodes and yield them per size
# Based on procedure RootedTree from https://dl.acm.org/doi/10.1145/1125994.1125995
def generate_trees(k_min, k_max):
    current_node = 1  # current node being generated
    current_trees = (
        []
    )  # contains <current_node> nodes trees generated after each iteration

    # first tree with only one node
    base = Node(current_node)
    current_trees.append(base)

    while current_node < k_max:
        if current_node >= k_min:
            yield current_trees

        # loop until every trees of <current_node> node were taken from the queue
        while current_trees[0].size() == current_node:
            tree = current_trees.pop(0)
            node = tree.get_node(current_node)

            # generate all possible trees of <current_node + 1> node from <tree>
            while node is not None:
                new_tree = tree.deep_copy()
                new_node = Node(current_node + 1, parent=new_tree.get_node(node.name))

                # check that the tree has not already been added
                if any(
                    isomorphic(new_tree, t)
                    for t in current_trees
                    if t.size() == current_node + 1
                ):
                    node = node.parent
                    continue

                current_trees.append(new_tree)
                node = node.parent

        current_node += 1

    yield current_trees


def label(parent, node, labels):
    grouped_children = []
    for c in node.children:
        for gc in grouped_children:
            if isomorphic(c, gc[0]):
                gc.append(c)
                break

    sizes = [sum(g.size() for g in gc) for gc in grouped_children]


def sized_partitions(data, sizes):
    assert len(data) == sum(sizes)

    if len(sizes) == 0:
        yield []

    for comb in combinations(data, sizes[0]):
        remaining = list(set(data).difference(set(comb)))
        for remaining_part in sized_partitions(
            remaining, sizes[1:] if len(sizes) > 1 else []
        ):
            yield [comb] + remaining_part


# for p in sized_partitions([1, 2, 3, 4], [1, 3]):
#     print(p)


"""
for i, batch in enumerate(get_trees(1, 5)):
    print("==================")
    print(f"TREES OF SIZE {i+1}")
    print("==================")
    for t in batch:
        print()
        t.print_tree()
"""

"""
root_1 = Node(10)
child1_1 = Node(11, parent=root_1)
child2_1 = Node(12, parent=root_1)
child3_1 = Node(13, parent=child1_1)
child4_1 = Node(14, parent=child1_1)

root_2 = Node(20)
child1_2 = Node(21, parent=root_2)
child2_2 = Node(22, parent=root_2)
child3_2 = Node(23, parent=child2_2)
child4_2 = Node(24, parent=child2_2)

root_3 = Node(30)
child1_3 = Node(31, parent=root_3)
child2_3 = Node(32, parent=child1_3)
child3_3 = Node(33, parent=root_3)
child4_3 = Node(34, parent=child3_3)

print("Tree 1 :")
root_1.print_tree()
print("Tree 2 :")
root_2.print_tree()
print("Tree 3 :")
root_3.print_tree()

print("T1 isomorphic to T2 ?", isomorphic(root_1, root_2))
print("T2 isomorphic to T1 ?", isomorphic(root_2, root_1))
print("T1 isomorphic to T3 ?", isomorphic(root_1, root_3))
"""
"""
# Create nodes
root = Node("Root")
level1_child1 = Node("level1_child1")
level1_child2 = Node("level1_child2")
level2_child1 = Node("level2_child1")
level2_child2 = Node("level2_child2")
level3_child1 = Node("level3_child1")
level3_child2 = Node("level3_child2")

# Attach level 1 nodes to root
level1_child1.attach(root)
level1_child2.attach(root)

# Attach level 2 nodes to level 1 nodes
level2_child1.attach(level1_child1)
level2_child2.attach(level1_child2)

# Attach level 3 nodes to level 2 nodes
level3_child1.attach(level2_child1)
level3_child2.attach(level2_child1)

# Now, print the entire tree
root.print_tree()
print()
level1_child1.print_tree()
print()
for n in [root, level1_child1, level1_child2, level2_child1, level2_child2, level3_child1, level3_child2]:
    print(n.name, [n.name for n in n.descendants])
print()
for n in [root, level1_child1, level1_child2, level2_child1, level2_child2, level3_child1, level3_child2]:
    print(n.name, n.depth)
print()

print("base tree")
root.print_tree()
print()

print("spr level2_child1 to root")
level2_child1.spr(root)
root.print_tree()
print()

print("spr level2_child1 to level2_child2")
level2_child1.spr(level2_child2)
root.print_tree()
print()

for n in [root, level1_child1, level1_child2, level2_child1, level2_child2, level3_child1, level3_child2]:
    print(n.name, [n.name for n in n.descendants])
print()
for n in [root, level1_child1, level1_child2, level2_child1, level2_child2, level3_child1, level3_child2]:
    print(n.name, n.depth)
print()

print("spr level2_child1 back to level1_child1")
level2_child1.spr(level1_child1)
root.print_tree()
print()

print("Copied tree")
root_copy = root.deep_copy()
root_copy.print_tree()
print()

level1_child1_copy = root_copy.get_node("level1_child1")
level1_child2_copy = root_copy.get_node("level1_child2")
level2_child1_copy = root_copy.get_node("level2_child1")
level2_child2_copy = root_copy.get_node("level2_child2")
level3_child1_copy = root_copy.get_node("level3_child1")
level3_child2_copy = root_copy.get_node("level3_child2")

print("spr level2_child1_copy to root_copy")
level2_child1_copy.spr(root_copy)
root_copy.print_tree()
print()

print()
for n in [root_copy, level1_child1_copy, level1_child2_copy, level2_child1_copy, level2_child2_copy, level3_child1_copy, level3_child2_copy]:
    print(n.name, [n.name for n in n.descendants])
print()
for n in [root_copy, level1_child1_copy, level1_child2_copy, level2_child1_copy, level2_child2_copy, level3_child1_copy, level3_child2_copy]:
    print(n.name, n.depth)

print("base tree")
root.print_tree()
print()

for n in [root, level1_child1, level1_child2, level2_child1, level2_child2, level3_child1, level3_child2]:
    print(n.name, [n.name for n in n.descendants])
print()
for n in [root, level1_child1, level1_child2, level2_child1, level2_child2, level3_child1, level3_child2]:
    print(n.name, n.depth)
print()

n1 = Node("appendedto1", parent=root.get_node("level2_child1"))
n2 = Node("appendedto2", parent=root.get_node("level3_child1"))
root.print_tree()
"""
