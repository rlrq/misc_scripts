import sys
sys.path.append("/mnt/chaelab/rachelle/src")
# from newick_extension import *

import basics

import newick
import time
import copy
import itertools


def set_leaf_attribute(t, a, f):
    """
    Set attribute (a) of leaf nodes as output of function (f) that
        takes the node as input.
    Arguments:
        t (newick.Node): tree read by package newick
        a (str): name of attribute in which to store the output of f
        f (function): function that takes a node and outputs a value
    Returns: None
    """
    nodes = t.get_leaves()
    for n in nodes:
        setattr(n, a, f(n))
    return

def set_leaf_attributes(t, a_f_dict):
    """
    Set multiple attributes for all leaf nodes in tree. Calls :func:`set_leaf_attribute`.
    Arguments:
        t (newick.Node): tree object read by package newick
        a_f_dict (dict): dictionary of {"<attribute name>": f}, where
            f is a function that takes a node and outputs a value to
            be stored in attribute with name <attribute name>.
    Returns: None
    """
    for a, f in a_f_dict.items():
        set_leaf_attribute(t, a, f)
    return

def filter_leaves(t, filter_function):
    """
    Apply filter function ('filter_function') and output a filtered list of nodes.
    """
    return [n for n in t.get_leaves() if filter_function(n)]

def prune_leaves_by_filter(t, filter_function, inverse = True,
                           remove_redundant_nodes = True, keep_leaf_name = True,
                           in_place = True):
    """
    Prune a tree based on filtered nodes.
        Calls :func:`filter_leaves_by_group` to obtain list of nodes filtered within
        an attribute group. See :func:`filter_leaves_by_group` for filtering details.
    Arguments:
        t (newick.Node): tree object read by package newick.
        filter_function (func): function that takes a list of nodes and
            outputs a filtered node/list of nodes
        inverse (bool): passed to t.prune(). Set to False to remove nodes output
            by :func:`filter_leaves_by_group`, True to keep. Default=True.
        remove_redundant_nodes (bool): collapse internal nodes with exactly 2
            vertices after pruning. Calls t.remove_redundant_nodes(). Default=True.
        keep_leaf_name (bool): passed to t.remove_redundant_nodes
        in_place (bool): modify input tree object directly. Set to False to create
            and output a deep copy. Default=True
    Returns:
        None if in_place=True
        tree of class newick.Node if in_place=False
    """
    ## create deep copy if requested
    if not in_place: t = copy.deepcopy(t)
    ## obtain list of filtered leaf nodes
    filtered_leaves = filter_leaves(t, filter_function)
    ## prune
    t.prune(filtered_leaves, inverse = inverse)
    ## collapse redundant nodes
    if remove_redundant_nodes:
        t.remove_redundant_nodes(keep_leaf_name = keep_leaf_name)
    return None if in_place else t

def filter_leaves_by_group(t, attribute_group, filter_function):
    """
    Group leaf nodes by an attribute (specified by 'attribute_group')
        and apply a filter function ('filter_function') that operates 
        within each attribute group and outputs a filtered node/list of
        nodes.
    Arguments:
        t (newick.Node): tree object read by package newick.
        attribute_group (str): name of an attribute which value will be
            used to group nodes for within-group filtering
        filter_function (func): function that takes a list of nodes and
            outputs a filtered node/list of nodes
    Returns:
        list of filtered nodes ([<newick.Node>])
    """
    ## group nodes by node.<attribute_group> value
    groups = basics.split_iter_by_key(t.get_leaves(), lambda n:getattr(n, attribute_group))
    ## filter within groups
    groups_filtered = {k: ([output]
                           ## standardise to iterable so we can use itertools.chain
                           if not (isinstance((output := filter_function(leaves)), list) or isinstance(output, tuple))
                           else output)
                       for k, leaves in groups.items()}
    ## ungroup and flatten to a single-level list
    return list(itertools.chain(*groups_filtered.values()))

def prune_leaves_by_group(t, attribute_group, filter_function, inverse = True,
                          remove_redundant_nodes = True, keep_leaf_name = True,
                          in_place = True):
    """
    Prune a tree based on filtered nodes.
        Calls :func:`filter_leaves_by_group` to obtain list of nodes filtered within
        an attribute group. See :func:`filter_leaves_by_group` for filtering details.
    Arguments:
        t (newick.Node): tree object read by package newick.
        attribute_group (str): name of an attribute which value will be
            used to group nodes for within-group filtering
        filter_function (func): function that takes a list of nodes and
            outputs a filtered node/list of nodes
        inverse (bool): passed to t.prune(). Set to False to remove nodes output
            by :func:`filter_leaves_by_group`, True to keep. Default=True.
        remove_redundant_nodes (bool): collapse internal nodes with exactly 2
            vertices after pruning. Calls t.remove_redundant_nodes(). Default=True.
        keep_leaf_name (bool): passed to t.remove_redundant_nodes
        in_place (bool): modify input tree object directly. Set to False to create
            and output a deep copy. Default=True
    Returns:
        None if in_place=True
        tree of class newick.Node if in_place=False
    """
    ## create deep copy if requested
    if not in_place: t = copy.deepcopy(t)
    ## obtain list of filtered leaf nodes
    filtered_leaves = filter_leaves_by_group(t, attribute_group, filter_function)
    ## prune
    t.prune(filtered_leaves, inverse = inverse)
    ## collapse redundant nodes
    if remove_redundant_nodes:
        t.remove_redundant_nodes(keep_leaf_name = keep_leaf_name)
    return None if in_place else t

def prune_unique_leaves(*trees, match_function = lambda n:n.name, in_place = True):
    """
    Only keep leaves that are shared by all input trees. Trees are modified in-place.
    Arguments:
        trees (newick.Node): tree object read by package newick.
        match_function (func): function by which to match leaves between trees
        in_place (bool): modify input tree objects directly. Set to False to create
            and output list of deep copies. Default=True
    Returns:
        None if in_place=True
        list of trees of class newick.Node if in_place=False
    """
    ## create deep copy if requested
    if not in_place:
        trees = [copy.deepcopy(t) for t in trees]
    ## per tree, group leaves by output of match_function
    match_values = [basics.invert_dict({n: match_function(n) for n in t.get_leaves()}) for t in trees]
    ## get match_function output(s) that are present in all trees
    valid_values = set.intersection(*map(lambda d:set(d.keys()), match_values))
    ## filter all trees to only keep leaves that have match_function outputs in all trees
    for i in range(len(match_values)):
        nodes_to_keep = list(itertools.chain(*[match_values[i].get(val, []) for val in valid_values]))
        trees[i].prune(nodes_to_keep, inverse = True)
    return None if in_place else trees

def is_single_vertex_tree(t):
    """
    'A tree is a single-vertex tree if it consists of a single vertex, which is a leaf of the tree.'
    Returns: bool
    """
    return t.ancestor is None and len(t.descendants) == 0

def is_single_edge_tree(t, disregard_root = True):
    """
    'A tree is a single-edge tree if it consists of a single edge, which contains two leaves and no non-leaf vertices.'
    Returns: bool
    """
    ## assume root is valid leaf; None -> t.ancestor (sole descendant is t) -> t -> None
    valid_1 = (not disregard_root and t.ancestor is not None and t.ancestor.ancestor is None and len(t.ancestor.descendants) == 1)
    ## assume root is valid leaf; None -> t -> t.descendant (N=1, is_leaf=True)
    valid_2 = (not disregard_root and t.ancestor is None and len(t.descendants) == 1 and t.descendants[0].is_leaf)
    ## disregard root node; None -> t (root) -> t.descendant (N=1, both is_leaf=True)
    valid_3 = (disregard_root and t.ancestor is None and len(t.descendants) == 2 and all(d.is_leaf for d in t.descendants))
    return valid_1 or valid_2 or valid_3

def sibling_set(n, disregard_root = True):
    """
    Gets list of sibling leaves (INcludes n) if input node n is a leaf AND it has sibling leaves.
    Else if n is an internal node, get its direct descendants that are leaves.
    'A sibling set is a set of leaves that are all siblings.'
    'Two leaves in a forest are siblings if they either share the same parent, or are the two leaves of a single-edge tree.'
    Arguments:
        n (newick.Node): node
    Returns: list of nodes (newick.Node)
    """
    ## if root is not counted as a leaf
    if disregard_root:
        if n.is_leaf: ## if n is leaf, get its siblings
            return [node for node in n.ancestor.descendants if node.is_leaf]
        else: ## if n is internal node, get its descendants that are leaves
            return [node for node in n.descendants if node.is_leaf]
    ## if root can be counted as leaf
    elif (not disregard_root and is_single_edge_tree(n)):
        return [node for node in {n}.union({n.ancestor}).union(set(n.descendants)) if node is not None]
    return []

def siblings(n, disregard_root = True):
    """
    Gets list of sibling leaves (EXcludes n) if input node n is a leaf AND it has sibling leaves.
    'Two leaves in a forest are siblings if they either share the same parent, or are the two leaves of a single-edge tree.'
    Arguments:
        n (newick.Node): node
    Returns: list of nodes (newick.Node)
    """
    return [node for node in sibling_set(n) if node is not n]

def is_sibling(n1, n2, disregard_root = True):
    """
    'Two leaves in a forest are siblings if they either share the same parent, or are the two leaves of a single-edge tree.'
    Returns: bool
    """
    shared_ancestor = n1.ancestor is n2.ancestor
    single_edge_tree = (is_single_edge_tree(n1, disregard_root = disregard_root) or
                        is_single_edge_tree(n2, disregard_root = disregard_root))
    return (n1 is not n2) and (shared_ancestor or single_edge_tree)

def is_bottommost_sibling_set(nodes, disregard_root = True):
    """
    'A bottommost sibling set (abbr. BSS) is a maximal sibling set X such that either the degree of their parent is at most
        |X|+1, or X is the leaf set of a single-edge tree.'
    """
    ## all nodes are leaves and the set of descendants of the first node's ancestor is equivalent to the set of nodes
    valid_1 = (all(map(lambda n:n.is_leaf, nodes)) and set(nodes[0].ancestor.descendants) == set(nodes))
    ## the nodes belong to a single-edge tree
    valid_2 = (not disregard_root and len(nodes) == 2 and is_single_edge_tree(node[0], disregard_root = False))
    valid_3 = (disregard_root and len(nodes) == 2 and is_single_edge_tree(node[0].ancestor, disregard_root = True))
    return valid_1 or valid_2 or valid_3

def get_root_node(n):
    """
    Get the root node of the tree in which n resides
    """
    if n.ancestor is not None:
        return root(n.ancestor)
    return n

def get_root(n):
    while n.ancestor is not None:
        n = n.ancestor
    return n

def match_across_forests(match_function = None):
    """Returns a function that takes a node and automatically selects a match function based on it"""
    if match_function is not None:
        return match_function
    def helper(n):
        ## get node's root node's forest
        forest = get_root(n).forest
        ## use the forest's resolve_match_function method
        return forest.resolve_match_function(n)
    return helper

class Forest():
    def __init__(self, trees = [], match_function = None):
        self.trees = [],
        self._match_function = match_function
    def match_function(self, t):
        """Get a function that matches leaves between trees. Uses tree's function if available, otherwise Forest's"""
        if t.match is not None:
            return t.match
        else:
            return self._match_function
    def resolve_match_function(self, n, match_function = None):
        """Get match function for a node. Uses node's tree's if available, otherwise Forest's"""
        if match_function is None:
            match_function = self.match_function(n.get_root_node())
        return match_function
    def match(self, n, match_function = None):
        """Function that takes a node and returns a value for matching between trees"""
        if match_function is None:
            match_function = self.match_function(n.get_root_node())
        return match_function(n)
    def filter_leaves(self, filter_function):
        """Returns list of nodes that return True when given to the filter_function"""
        output = []
        for t in self.trees:
            output.extend(t.filter_leaves(filter_function = filter_function))
        return output
    def drop_leaf_by_match(self, value, match_function = None, keep_leaf = True, inverse = False, **kwargs):
        """Takes a value and drops/separates leaves that return that value when given to a match function"""
        new_trees = []
        drop_trees = []
        for i, t in self.trees:
            match_tree = t.match
            match = self.resolve_match_function(t, match_function)
            ## leaves to drop
            matched_leaves = t.filter_leaves(filter_function = (lambda n:(match(n)==value)!=inverse))
            if not matched_leaves: continue
            ## handle single-vertex trees
            if t.is_single_vertex_tree():
                if not keep_leaf: drop_trees.append(i)
                continue
            new_trees.extend(matched_leaves)
            ## if not dropping leaves from forest, we need to update new trees' roots with metadata
            ## default behaviour; cleaved leaves become new, single-vertex trees
            if keep_leaf:
                ## populate forest + match function info for new trees
                for new_t in new_trees:
                    new_t.match = match_tree
                    new_t.forest = self
            ## filtering is in-place
            matched_leaves = set(matched_leaves)
            t.prune_leaves_by_filter(filter_function = (lambda n:n in matched_leaves),
                                     in_place = True, inverse = False, keep_leaf_name = True,
                                     **{k:v for k, v in kwargs.items() if k not in {"in_place", "keep_leaf_name"}})
        ## default behaviour; cleaved leaves become new, single-vertex trees
        if keep_leaf:
            self.trees.extend(new_trees)
        else: ## alternative behaviour, previously single-vertex trees with dropped leaves are removed from Forest
            ## tbh drop_trees would already be empty if keep_leaf=False but we're still gonna check anyway
            for i in drop_trees[::-1]: ## reverse order for popping so we pop from end; avoids messing up indexing
                _ = self.trees.pop(i)
        return
    def reduction_1(self, f2, match_function = None, **kwargs):
        """
        Reduction Rule 1: If a label l is in a single-vertex tree in one of the
        forests F1 or F2, then remove the edge (if any) incident to the label l
        in the other forest.
        """
        _match_function = match_function
        for t in self.trees:
            if match_function is None:
                if t.match_function is not None: _match_function = t.match
                else: _match_function = self.match
            if t.is_single_vertex_tree():
                match_value = _match_function(t)
                f2.drop_leaf_by_match(match_value, match_function = match_function, keep_leaf = True,
                                      remove_redundant_nodes = True, **kwargs)
        ## for proper reduction, this function needs to also be called by f2
        ## (ideally keep calling by self & f2 alternately until number of trees in both forests remains unchanged
        return
    def reduction_1_reciprocal(self, f2, **kwargs):
        """Performs reduction_1 on both trees alternately until number of trees stabilises"""
        num_trees_1 = len(self.trees)
        num_trees_2 = len(f2.trees)
        while True:
            self.reduction_1(f2, **kwargs)
            f2.reduction_1(self, **kwargs)
            new_num_trees_1 = len(self.trees)
            new_num_trees_2 = len(f2.trees)
            ## if reduction no longer changes number of trees
            if num_trees_1 == new_num_trees_1 and num_trees_2 == new_num_trees_1:
                break
            else: ## else update number of trees
                num_trees_1 = new_num_trees_1
                num_trees_2 = new_num_trees_2
        return
    def shrink(self, nodes):
        """
          By shrinking X, we mean we delete all leaves in X, and introduce a new leaf
        v_x with a new label l_x (e.g. we can use l(X) for l_x) and let v_x be adjacent
        to the common neighbour of th leaves in X if such a common neighbour exists.
        """
        pass
    def reduction_2(self, f2, x1, match_function = None, disregard_root = True, **kwargs):
        """
          Let X_1 be a BSS in F_1, and let X_2 be a set of leaves in F_2 such that
        l(X_1) = l(X_2). If X_2 is the leaf set of a single-edge tree in F_2, or if
        X_2 is a sibling set whose parent has degree at most |X_2|+1, then shrink
        X_1 in F_1, and shrink X_2 in F_2.
          Note that after applying Reduction Rule 2, if a vertex v in either F_1 or
        F_2 is adjacent to the new leaf (there is at most one such vertex), then the
        degree of v is at most 2. In particular, if v is unlabelled, then v will be
        contracted.
        """
        ## check if x1 is BSS
        if not is_bottommost_sibling_set(x1, disregard_root = disregard_root):
            print("X1 is not a BSS. Exiting function.")
            return
        ## get x1 values for matching
        x1_values = set(match_across_forests(match_function = match_function)(n)(n) for n in x1)
        ## get equivalent nodes in f2
        x2 = f2.filter_leaves(filter_function = lambda n:match_across_forests(match_function = match_function)(n)(n))
        ## if x2 is also a BSS, collapse in both trees
        if is_bottommost_sibling_set(x2, disregard_root = disregard_root):
            ## shrink x1 and x2 leaves into new node l
            l1 = ShrunkenTree(x1)
            l2 = ShrunkenTree(x2)
            ## create common name for these new nodes
            match_value = str(id(l1)) + '_' + str(id(l2))
            l1._name = match_value
            l2._name = match_value
            ## store each other as equivalents
            l1.equivalents.append(l2)
            l2.equivalents.append(l1)
            ## get common neighbour of the leaves in each tree
            ancestor1 = x1[0].ancestor
            ancestor2 = x2[0].ancestor
            ## remove X1 an X2 from original tree
            ancestor1.prune(x1)
            ancestor2.prune(x2)
            ## replace X1 and X2 qith l1 and l2
            ancestor1.add_descendant(l1)
            ancestor2.add_descendant(l2)
        return
    def is_subgraph(self, f2):
        """Returns True if self is subgraph of f2"""
        ## use graph_tool intead
        pass

class ShrunkenTree(newick.Node):
    def __init__(self, nodes, equivalents = [], **kwargs):
        super().__init__()
        self.nodes = list(nodes) ## nodes can also include ShrunkenTree objects
        self.equivalents = equivalents


## update class with new methods
newick.Node.set_leaf_attribute = set_leaf_attribute
newick.Node.set_leaf_attributes = set_leaf_attributes
newick.Node.filter_leaves_by_group = filter_leaves_by_group
newick.Node.prune_leaves_by_group = prune_leaves_by_group
newick.Node.filter_leaves = filter_leaves
newick.Node.is_single_vertex_tree = is_single_vertex_tree
newick.Node.is_single_edge_tree = is_single_edge_tree
newick.Node.sibling_set = sibling_set
newick.Node.siblings = siblings
newick.Node.is_sibling = is_sibling
newick.Node.root = get_root_node
newick.Node.forest = None
newick.Node.match = None


def prep_for_edgelist(t):
    v = []
    e = []
    def recur(node):
        for child in node.descendants:
            v.append(child)
            e.append([node, child])
            recur(child)
        return
    v.append(t)
    recur(t)
    return v, e

def string_to_double(s):
    try:
        s = float(s)
    except ValueError:
        return ''
    except TypeError:
        return ''
    return s

default_vertex_property_map = {
    "bootstrap": string_to_double
}
default_edge_property_map = {
    "length": lambda src, dst: dst.length
}
default_vertex_label_function = lambda v:v.name
default_edge_label_function = lambda src,dst:''

def format_string_v1(v):
    return v if not isinstance(v, str) else f'"{v}"'

def write_edgelist(fout, v_list, e_list,
                   v_label_function = default_vertex_label_function,
                   e_label_function = default_edge_label_function,
                   v_property_map = default_vertex_property_map,
                   e_property_map = default_edge_property_map,
                   v_property_order = None, e_property_order = None):
    v_dict = {}
    v_property_order = v_property_order if v_property_order is not None else sorted(v_property_map.keys())
    e_property_order = e_property_order if e_property_order is not None else sorted(e_property_map.keys())
    ## string values should be in double quotes
    v_property_order_reformat = [format_string_v1(name) for name in v_property_order]
    e_property_order_reformat = [format_string_v1(name) for name in e_property_order]
    with open(fout, "w+") as f:
        ## write property names in order
        _ = f.write(f"# vertex {' '.join(v_property_order_reformat)}\n")
        _ = f.write(f"# edge {' '.join(e_property_order_reformat)}\n")
        ## write all vertices
        for i, v in enumerate(v_list):
            v_dict[v] = i
            ## string values should be in double quotes
            v_label = format_string_v1(v_label_function(v))
            property_dict = {name: format_string_v1(f(v)) for name, f in v_property_map.items()}
            _ = f.write(f"{i} * {v_label} {' '.join(str(property_dict[name]) for name in v_property_order)}\n")
        ## write all edges
        for src, dst in e_list:
            ## string values should be in double quotes
            e_label = format_string_v1(e_label_function(src, dst))
            property_dict = {name: format_string_v1(f(src, dst)) for name, f in e_property_map.items()}
            _ = f.write(f"{v_dict[src]} {v_dict[dst]} {e_label} {' '.join(str(property_dict[name]) for name in e_property_order)}\n")
    return

def newick_to_edgelist(fin, fout, **kwargs):
    print("---Reading newick file---")
    t = newick.read(fin)[0]
    print("---Parsing vertices and edges---")
    v_list, e_list = prep_for_edgelist(t)
    print("---Writing edgelist file---")
    write_edgelist(fout, v_list, e_list, **kwargs)
    return

def parse_edgelist(fin, v_property_coerce = {}, e_property_coerce = {}):
    import csv
    v_properties = []
    e_properties = []
    v_dict = {}
    e_list = []
    v_coerce = []
    e_coerce = []
    with open(fin, 'r') as f:
        reader = csv.reader(f, delimiter = ' ', quotechar = '"')
        for line in reader:
            ## if vertex
            if line[1] == '*':
                ## index vertex by unique integer identifier
                v_dict[int(line[0])] = [v_coerce[i](x) for i, x in enumerate(line[2:])]
            ## if edge
            elif line[0] != '#':
                ## index edges by unique integer identifiers for src and dst
                e_list.append([int(line[0]), int(line[1])] + [e_coerce[i](x) for i, x in enumerate(line[2:])])
            ## if comment
            else:
                if line[1] == "vertex":
                    v_properties = ["label"] + line[2:]
                    v_coerce = [v_property_coerce.get(name, lambda x:x) for name in v_properties]
                elif line[1] == "edge":
                    e_properties = ["label"] + line[2:]
                    e_coerce = [e_property_coerce.get(name, lambda x:x) for name in e_properties]
    return v_properties, e_properties, v_dict, e_list

def edgelist_to_graph(fin, **kwargs):
    import graph_tool.all as gt
    v_properties, e_properties, v_dict, e_list = parse_edgelist(fin, **kwargs)
    eprops = [(eprop, reformat_type(infer_type(e[i+2] for e in e_list).__name__))
              for i, eprop in enumerate(e_properties)]
    g = gt.Graph(e_list, eprops = eprops)
    ## assign vertex properties
    for vprop_i, v_prop in enumerate(v_properties):
        ## creater property
        vprop = g.new_vertex_property(reformat_type(infer_type(v[vprop_i] for v in v_dict.values()).__name__))
        g.vp[v_prop] = vprop
        ## iterate through vertices and assign property values
        for v_i, v in v_dict.items():
            vert = g.vertex(v_i) ## get vertex
            vprop[vert] = v[vprop_i] ## assign
    ## add collapsed_vertices
    vprop = g.new_vertex_property("object")
    g.vp["collapsed_vertices"] = vprop
    for v_i in v_dict:
        vert = g.vertex(v_i)
        vprop[vert] = []
    return g

def infer_type(l):
    type_count = basics.get_count_dict(map(type, l))
    valid_types = [t for t in type_count if t is not type(None)]
    max_type = max(type_count, key = lambda k:type_count[k])
    return max_type

def reformat_type(v):
    if v == "str": return "string"
    elif v == "float": return "double"
    elif v in ["dict", "list"]: return "object"
    return v
