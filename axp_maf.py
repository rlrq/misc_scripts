import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import newick

## based on https://www.sciencedirect.com/science/article/pii/S0304397514008238#br0060
## Parameterized and approximation algorithms for maximum agreement forest in multifurcating trees
## 2015; Jianer Chen, Jia-Hao Fan, & Sing-Hoi Sze; Theoretical Computer Science. 562:496-512

#### for development, use these trees, which have some multifurcating subtrees
from Bio import SeqIO
import copy
import itertools
import regex as re
import basics
## ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1.nwk has an accession with 2 genes; one is nested within the other, keep the longer one
nwk1 = "/mnt/chaelab/rachelle/nlr_survey/results/RNL/tree/ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1.nwk"
## ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1-L2.nwk has genes with multiple splice variants; needs to be processed to select the longest variant per accession
nwk2 = "/mnt/chaelab/rachelle/nlr_survey/results/RNL/tree/ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1-L2.nwk"
## to select longest sequence, we need to know sequence lengths, so read in the alignment, drop gaps, get length
aln = "/mnt/chaelab/rachelle/nlr_survey/results/RNL/aln/ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.fasta"
seqlens = {rec.id: len(rec.seq.replace('-', '')) for rec in SeqIO.parse(aln, "fasta")}

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

def is_bottommost_sibling_set(*nodes, disregard_root = True):
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

class Forest():
    def __init__(self, trees = [], match_function = None):
        self.trees = [],
        self._match_function = match_function
    def match_function(self, t):
        if t.match is not None:
            return t.match
        else:
            return self._match_function
    def resolve_match_function(self, n, match_function = None):
        if match_function is None:
            match_function = self.match_function(n.get_root_node())
        return match_function
    def match(self, n, match_function = None):
        if match_function is None:
            match_function = self.match_function(n.get_root_node())
        return match_function(n)
    def drop_leaf_by_match(self, value, match_function = None):
        for t in self.trees:
            match = self.resolve_match_function(t, match_function)
            t.prune_leaves_by_filter(match_function = match)
            pass
        pass
    def reduction_1(self, f2, match_function = None):
        """
        Reduction Rule 1: If a label l is in a single-vertex tree in one of the
        forests F1 or F2, then remove the edge (if any) incident to the label l
        in the other forest.
        """
        for t in self.trees:
            if match_function is None:
                if t.match_function is not None: match_function = t.match
                else match_function = self.match
            if t.is_single_vertex_tree():
                match_value = match_function(t)
                

## update class with new methods
newick.Node.set_leaf_attribute = set_leaf_attribute
newick.Node.set_leaf_attributes = set_leaf_attributes
newick.Node.filter_leaves_by_group = filter_leaves_by_group
newick.Node.prune_leaves_by_group = prune_leaves_by_group
newick.Node.is_single_vertex_tree = is_single_vertex_tree
newick.Node.is_single_edge_tree = is_single_edge_tree
newick.Node.sibling_set = sibling_set
newick.Node.siblings = siblings
newick.Node.is_sibling = is_sibling
newick.Node.root = get_root_node

## some functions for assigning attributes + filtering
get_gid = lambda seqid: re.search("(?<=\\|).+(?=\\|CDS\\|)", seqid).group(0)
get_accid = lambda seqid: ("TAIR10" if "Col-0_ref" in seqid
                           else re.search("(?<=Reference\\|)[^|]+", seqid).group(0) if re.search("\\d\\|G\\d", seqid)
                           else re.search("(?<=Reference\\|).+(?=\\.g[-\\d]+?\\|CDS)", seqid).group(0))
attr_map = {"seqlen": lambda n:seqlens[n.name],
            "accid": lambda n:get_accid(n.name),
            "gid": lambda n:get_gid(n.name)}

## read & process trees
t1 = newick.read(nwk1)[0]
t1.set_leaf_attributes(attr_map)
t1.prune_leaves_by_group("accid", lambda leaves:max(leaves, key = lambda leaf:leaf.seqlen), in_place = True)
t1.match = lambda n:n.accid
t2 = newick.read(nwk2)[0]
t2.set_leaf_attributes(attr_map)
t2.prune_leaves_by_group("accid", lambda leaves:max(leaves, key = lambda leaf:leaf.seqlen), in_place = True)
prune_unique_leaves(t1, t2, match_function = lambda n:n.accid, in_place = True)
t2.match = lambda n:n.accid


###################
##  DEFINITIONS  ##
###################

## (direct quotes from paper) ##

# In this paper, all graphs are undirected.
# All trees in our discussion are unrooted.

## ---------------------------------------------------------------------------
## --* leaf
## --* forest
## --** v e E' G F L l 
## ---------------------------------------------------------------------------
#   For a vertex v, an edge e, and an edge subset D' in a graph G, denote by
# G-v, G-e, and G-E' the graphs obtained from G with v, e, and the edges E'
# removed, respectively. A leaf of a tree is a vertex of degree less than 2.
# A forest is a collection of disjoint trees. A nonempty forest F is leaf-
# labelled over a label-set L if there is a one-to-one mapping from the leaves
# of F to the elements of L (with all non-leaf vertices unlabelled). The label
# of a leaf v is denoted by l(v). More generally, for a subforest F' of F,
# denote by l(F') the set of labels for the leaves in F'.

## ---------------------------------------------------------------------------
## --* single-vertex tree
## --* single-edge tree
## ---------------------------------------------------------------------------
#   A tree is a single-vertex tree if it consists of a single vertex, which is
# a leaf of the tree. A tree is a single-edge tree if it consists of a single
# edge, which contains two leaves and no non-leaf vertices. The parent of a
# leaf v in a tree with at least three vertices is the unique non-leaf vertex
# adjacent to v.

## ---------------------------------------------------------------------------
## --* siblings
## --* bottommost sibling set (BSS)
## ---------------------------------------------------------------------------
#   Two leaves in a forest are siblings if they either share the same parent,
# or are the two leaves of a single-edge tree. A sibling set is a set of
# leaves that are all siblings. A bottommost sibling set (abbr. BSS) is a
# maximal sibling set X such that either the degree of their parent is at most
# |X|+1, or X is the leaf set of a single-edge tree. By definition, the leaf
# of a single-vertex tree is not a BSS. Moreover, by our assumption, an
# unlabelled vertex has degree at least 3. Thus, a BSS contains at least two
# leaves, and a leaf-labelled tree that is not a single-vertex must contain a
# BSS.

def reduction_rule_1(t, disregard_root = True):
    """
      If a label l is in a single-vertex tree in one of the forests F_1 and F_2,
    then remove the edge (if any) incident to the label l in the other forest.
    """
    pass

def reduction_rule_2(t):
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
    pass

