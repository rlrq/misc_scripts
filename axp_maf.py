import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from newick_extension import *

# import newick
# import time
# import copy

## based on https://www.sciencedirect.com/science/article/pii/S0304397514008238#br0060
## Parameterized and approximation algorithms for maximum agreement forest in multifurcating trees
## 2015; Jianer Chen, Jia-Hao Fan, & Sing-Hoi Sze; Theoretical Computer Science. 562:496-512

#### for development, use these trees, which have some multifurcating subtrees
from Bio import SeqIO
# import itertools
import regex as re
# import basics
## ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1.nwk has an accession with 2 genes; one is nested within the other, keep the longer one
nwk1 = "/mnt/chaelab/rachelle/nlr_survey/results/RNL/tree/ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1.nwk"
## ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1-L2.nwk has genes with multiple splice variants; needs to be processed to select the longest variant per accession
nwk2 = "/mnt/chaelab/rachelle/nlr_survey/results/RNL/tree/ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.rooted.subtree-ADR1-L2.nwk"
## to select longest sequence, we need to know sequence lengths, so read in the alignment, drop gaps, get length
aln = "/mnt/chaelab/rachelle/nlr_survey/results/RNL/aln/ADR1-NRG1.col0-vdw-wlodzimierz.curated.pep.mafft.fasta"
seqlens = {rec.id: len(rec.seq.replace('-', '')) for rec in SeqIO.parse(aln, "fasta")}

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

## convert to format that graph_tool can read
import networkx
import io
from Bio import Phylo
def newick_to_graphml(newick_string, fout):
    phylo = Phylo.read(io.StringIO(newick_string), "newick")
    g = Phylo.to_networkx(phylo)
    networkx.write_graphml(g, fout)
    return

t1_snwk = newick.dumps(t1)
t2_snwk = newick.dumps(t2)
graphml_1 = "/mnt/chaelab/rachelle/tmp/t1.graphml"
graphml_2 = "/mnt/chaelab/rachelle/tmp/t2.graphml"
## bootstrap values are converted to weight, I believe
## branch length info is lost
newick_to_graphml(t1_snwk, graphml_1)
newick_to_graphml(t2_snwk, graphml_2)

## convert to edge list (text file, format based on https://docs.oracle.com/en/database/oracle/property-graph/23.1/spgdg/edge-list-edge_list.html)
def parse_if_non_numeric(v, f):
    try:
        v = float(v)
        return ''
    except ValueError:
        return f(v)
    except TypeError:
        return ''

edgelist_1 = "/mnt/chaelab/rachelle/tmp/t1.edgelist"
edgelist_2 = "/mnt/chaelab/rachelle/tmp/t2.edgelist"
v_attr_map = {
    "seqlen": lambda n:seqlens.get(n.name, ''),
    "accid": lambda n:parse_if_non_numeric(n.name, get_accid),
    "gid": lambda n:parse_if_non_numeric(n.name, get_gid),
    "leaf": lambda n:n.is_leaf,
    "collapsed_vertices": lambda n:[]}
edgelist_kwargs = {
    "v_label_function": lambda v:v.name,
    "e_label_function": lambda src,dst:'',
    "v_property_map": v_attr_map,
    "e_property_map": {"length": lambda src,dst:dst.length}}
v_property_coerce = {
    "seqlen": lambda x:(-1 if x == '' else int(x)),
    "accid": lambda x:(None if x == '' else x),
    "gid": lambda x:(None if x == '' else x),
    "leaf": lambda x:(True if (low:=x.lower())=="true" else False if low=="false" else None)}
e_property_coerce = {
    "length": lambda x:(float("NaN") if x == '' else float(x))}
vcoll_attr_map = {
    "seqlen": lambda g,v:-1,
    "accid": lambda g,v:(';'.join(sorted(get_vprop_recursively(g, v, "accid")))),
    "gid": lambda g,v:(';'.join(sorted(get_vprop_recursively(g, v, "gid")))),
    "leaf": lambda g,v:True,
    "collapsed_vertices": lambda g,n:[]
}
ecoll_attr_map = {
    "length": lambda g,e:-1
}
econt_attr_map = {
    "length": lambda g,v:sum(g.ep["length"][e] for e in v.all_edges())
}
## convert to edgelist
tmp1 = "/mnt/chaelab/rachelle/tmp1.txt"
tmp2 = "/mnt/chaelab/rachelle/tmp2.nwk"
newick.write(t1, tmp1)
newick.write(t2, tmp2)
newick_to_edgelist(tmp1, edgelist_1, **edgelist_kwargs)
newick_to_edgelist(tmp2, edgelist_2, **edgelist_kwargs)
## read edgelist into graph-tool graph
g1 = edgelist_to_graph(edgelist_1, v_property_coerce = v_property_coerce, e_property_coerce = e_property_coerce)
g2 = edgelist_to_graph(edgelist_2, v_property_coerce = v_property_coerce, e_property_coerce = e_property_coerce)
## undirected
g1.set_directed(False)
g2.set_directed(False)

## update gt.Vertex object with new method
import graph_tool.all as gt
gt.Vertex.total_degree = lambda v:(v.in_degree() + v.out_degree())
# gt.Vertex.collapsed_vertices = []
gt.Vertex.summarise_collapsed_vertices = lambda v,g,vprop:';'.join(str(g.vp[vprop][vcoll]) for vcoll in g.vp["collapsed_vertices"][v])

def get_vprop_recursively(g, vcoll, attr):
    if not g.vp["collapsed_vertices"][vcoll]:
        return [g.vp[attr][vcoll]]
    output = []
    for v in g.vp["collapsed_vertices"][vcoll]:
        output.append(get_vprop_recursively(v, attr))
    return output

def is_single_vertex_tree(v):
    """A tree is a single-vertex tree if it consists of a single vertex, which is a leaf of the tree."""
    return v.total_degree() == 0

def is_single_edge_tree(v_e):
    """A tree is a single-edge tree if it consists of a single edge, which contains two leaves and no non-leaf vertices."""
    if isinstance(v_e, gt.Vertex):
        neighbours = list(v_e.all_neighbours())
        one_neighbour = len(neighbours) == 1
        if not one_neighbour: return False
        return neighbours[0].total_degree() == 1
    elif isinstance(v_e, gt.Edge):
        return v_e.source().total_degree() == v_e.target().total_degree() == 1
    raise Exception("Input must be a graph_tool vertex or edge object.")

def in_same_subgraph(g, v_list):
    """Returns true if vertices in list are in same connected subgraph (i.e. component)"""
    vertex_indices = [g.vertex_index[v] for v in v_list]
    g_is_directed = g.is_directed()
    g.set_directed(False)
    components = gt.label_components(g)[0].a
    g.set_directed(g_is_directed)
    return len(set(components[i] for i in vertex_indices)) == 1

def vertices_of_single_edge_tree(g, v_list):
    """Returns true if vertices are from a single-edge tree. Assumes single-edge tree has exactly 1 edge and 2 vertices."""
    if not in_same_subgraph(g, v_list): return False
    if len(set(v_list)) != 2: return False
    neighbours = set(itertools.chain(*[v.all_neighbours() for v in v_list]))
    return neighbours == set(v_list)

def reduction_rule_1(F1, F2, vp_compare):
    """
      If a label l is in a single-vertex tree in one of the forests F_1 and F_2,
    then remove the edge (if any) incident to the label l in the other forest.
    """
    ## extract vp_compare values from vertices in single vertex trees
    def single_vertex_attr(f):
        single_vertex_trees = []
        for v in f.vertices():
            if is_single_vertex_tree(v) and f.vp["leaf"][v]:
                single_vertex_trees.append(f.vp[vp_compare][v])
        return single_vertex_trees
    ## delete all edges adjacet to vertices whose vp_compare value is in vprop_val
    def prune(f, vprop_val):
        vprop_val = set(vprop_val)
        edges_removed = 0
        for v in f.vertices():
            if f.vp[vp_compare][v] in vprop_val:
                edges = v.all_edges()
                for e in edges:
                    f.remove_edge(e)
                    edges_removed += 1
        return edges_removed ## return this so we can determine when to stop
    ## reduce repeatedly until no more reduction is possible
    while True:
        F1_single_vertex_prop = single_vertex_attr(F1)
        F2_single_vertex_prop = single_vertex_attr(F2)
        p1 = prune(F1, F2_single_vertex_prop)
        p2 = prune(F2, F1_single_vertex_prop)
        ## exit if no reduction occurs in this iteration
        if p1 == p2 == 0: break
    return

def reduction_rule_2(F1, F2, X1, vp_compare,
                     vcoll_property_map = {}, ecoll_property_map = {},
                     econt_property_map = {}):
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
    if not is_bss(F1, X1): raise Exception("X1 is not a BSS.")
    ## get equivalent X2 vertices
    vp_X1 = set(F1.vp[vp_compare][v] for v in X1)
    X2 = [v for v in F2.vertices() if F2.vp[vp_compare][v] in vp_X1]
    if not is_bss(F2, X2): return ## nothing to reduce if X2 is not a BSS
    ## shrink
    vcoll1 = list(shrink(F1, X1, v_property_map = vcoll_property_map, e_property_map = ecoll_property_map).all_neighbours())
    vcoll2 = list(shrink(F2, X2, v_property_map = vcoll_property_map, e_property_map = ecoll_property_map).all_neighbours())
    ## contract (if necessary)
    v1 = None if len(v:=list(vcoll1.all_neighbours())) == 0 else v[0]
    v2 = None if len(v:=list(vcoll2.all_neighbours())) == 0 else v[0]
    if v1: _ = contract_vertex(F1, v1, e_property_map = econt_property_map)
    if v2: _ = contract_vertex(F2, v2, e_property_map = econt_property_map)
    # ## contract repeatedly until no further contractions are possible
    # while True: ## contract returns the number of vertices that were contracted
    #     c1 = contract(F1)
    #     c2 = contract(F2)
    #     if c1 == c2 == 0: break
    return

def shrink(g, x, v_property_map = {}, e_property_map = {}):
    """Replace vertices in X with a single vertex BUT ONLY IF X is a BSS."""
    if not is_bss(g, x): return
    parent = get_parent(g, x)
    ## create new vertex
    new_v = g.add_vertex()
    g.vp["collapsed_vertices"][new_v] = list(set(x)) ## store X as collapsed_vertices
    ## set vertex properties
    for name, f in v_property_map.items():
        g.vp[name][new_v] = f(g, new_v)
    ## connect vertex to parent (if parent exists)
    if parent:
        new_e = g.add_edge(parent, new_v)
        ## set edge properties
        for name, f in e_property_map.items():
            g.ep[name][new_e] = f(g, new_e)
    ## remove edges between parent and vertices in X
    for v in x:
        e = e1 if (e1:=g.edge(parent, v)) is not None else g.edge(v, parent)
        g.remove_edge(e)
    return new_v ## return the new vertex

def contract_vertex(g, v, e_property_map = {}):
    """
    Drop internal node with total degree 2, connect its neighbours with a new edge. For undirected graphs only.
    Return number of contractions.
    """
    ## skip if leaf
    if g.vp["leaf"][v]:
        return 0
    ## skip if total degree != 2
    elif v.total_degree() != 2:
        return 0
    else:
        neighbours = list(v.all_neighbours())
        ## create new edge between v's neighbours
        new_e = g.add_edge(*neighbours)
        ## set edge properties
        for name, f in e_property_map.items():
            g.ep[name][new_e] = f(g, v)
        ## delete vertex
        g.remove_vertex(v)
        return 1

def contract(g, e_property_map = {}):
    """Drop internal nodes with total degree 2, connect its neighbours with a new edge. Works properly for undirected graphs only."""
    num_contractions = 0
    ## vertex removal updates index of every vertex after it
    ## so we shouldn't remove while iterating
    ## and for efficiency we should also remove from decreasing index order
    v_list = reversed(sorted(g.vertices()))
    for v in v_list:
        num_contractions += contract_vertex(g, v, e_property_map = e_property_map)
    return num_contractions ## return the number of contractions so that we can determine when to stop contracting

def get_parent(g, x):
    if not is_sibling_set(g, x): raise Exception("X is not a sibling set.")
    parents = list(set.intersection(*[set(v.all_neighbours()) for v in x]))
    if len(parents) != 1:
        print("Warning: X has more/less than one common parent.")
        return parents
    return parents[0]

def is_bss(g, x):
    """
    'A bottommost sibling set (abbr. BSS) is a maximal sibling set X such that either the degree of their parent is at most
        |X|+1, or X is the leaf set of a single-edge tree.'
    """
    if vertices_of_single_edge_tree(g, x): return True ## return True if x's leaves are from the same single-edge tree
    if not in_same_subgraph(g, x): return False ## return False if vertices are from disjointed trees
    x = set(x)
    all_leaves = all(map(lambda v:(v.total_degree() == 1), x))
    if not all_leaves: return False
    neighbours = set(itertools.chain(*[v.all_neighbours() for v in x]))
    if len(neighbours) != 1: return False ## return False if leaves don't all share the same parent
    return list(neighbours)[0].total_degree() <= (len(x)+1) ## return whether parent is at most |X|+1

def is_sibling_set(g, x):
    """
    'A sibling set is a set of leaves that are all siblings.'
    'Two leaves in a forest are siblings if they either share the same parent, or are the two leaves of a single-edge tree.'
    """
    if vertices_of_single_edge_tree(g, x): return True ## return True if x's leaves are from the same single-edge tree
    if not in_same_subgraph(g, x): return False ## return False if vertices are from disjointed trees
    x = set(x)
    all_leaves = all(map(lambda v:(v.total_degree() == 1), x))
    if not all_leaves: return False
    neighbours = set(itertools.chain(*[v.all_neighbours() for v in x]))
    return len(neighbours) == 1 ## return whether siblings share the same parent

def vertices(g, *indices):
    return [g.vertex(i) for i in indices]

def assign_leaf(g):
    leaf = g.new_vertex_property("bool")
    for v in g.vertices():
        leaf[v] = v.out_degree() == 0 if g.is_directed else v.total_degree() == 1
    g.vp["leaf"] = leaf
    return

## some tests
g = gt.Graph([(1,2),(1,3),(1,9),(0,1),(0,4),(4,5),(4,6),(9,7),(7,8),(7,10)])
assign_leaf(g)

## randomly assign edge length
import random
l = g.new_edge_property("int")
g.ep["length"] = l
for e in g.edges():
    l[e] = random.randint(1,5)

cv = g.new_vertex_property("object")
g.vp["collapsed_vertices"] = cv
for v in g.vertices():
    cv[v] = []

ss = vertices(g, 2,3) ## sibling set (not bss)
bss = vertices(g, 2,3,9) ## bss
not_ss = vertices(g, 1,2,3)
dj = vertices(g, 0,7) ## disjointed

is_single_edge_tree(g.vertex(7)) ## True
is_sibling_set(g, ss)
is_bss(g, bss)

tmp = g.vertex(10) ## last vertex
g.remove_vertex(g.vertex(9))
tmp ## <Vertex object with index '10' at 0x7fdcfd94d540>

tmp = g.vertex(5)
tmp ## <Vertex object with index '5' at 0x7fdcfd94d540>
g.remove_vertex(g.vertex(4))
tmp ## <Vertex object with index '5' at 0x7fdcfd94d540>
g.vertex(4) ## <Vertex object with index '4' at 0x7fdcfd94d7c0> (vertex at index 4 should be the vertex that was in position 5)
## seems like the vertex object doesn't store properties of a given vertex, but just the index. when indices are updated, the object just points to whatever vertex is at that given index now
## we can't rely on them to keep track of the same vertex after re-indexing (unless we try to hash it?)
## honestly it even seems like every time we call g.vertex(n) we get a different object even with the same n and without any re-indexing

# g = gt.Graph([('1','2'),('1','3'),('1','9'),('0','1'),('0','4'),('4','5'),('4','6'),('9','7'),('7','8'),('7','10')], hashed = True)
g = gt.Graph([('6','0'),('0','1'),('1','2'),('2','3'),('3','4'),('4','5'),('5','6')], hashed = True)
ogi = g.new_vertex_property("int")
g.vp["ogi"] = ogi
for i,v in enumerate(g.vertices()):
    ogi[v] = i

tmp = g.vertex(6) ## last vertex?
g.vp["ogi"][tmp] ## 6
tmp2 = g.vertex('6')
g.vp["ogi"][tmp2] ## 6
g.vp.ids[tmp2]
g.remove_vertex(5)
tmp ## invalid vertex
g.vertex_index[tmp]
tmp = g.vertex(7) ## <Vertex object with index '7' at 0x7fdcfd94d7c0>
tmp2 = g.vertex('7') ## <Vertex object with index '7' at 0x7fdcfd94d640>
g.remove_vertex(6)
tmp ## <Vertex object with index '7' at 0x7fdcfd94d7c0>
tmp2 ## <Vertex object with index '7' at 0x7fdcfd94d640>
## sigh even if we retrieve by index 


g = gt.Graph([('six','zero'),('zero','one'),('one','two'),('two','three'),('three','four'),('four','five'),('five','six')], hashed = True)
ogi = g.new_vertex_property("int")
g.vp["ogi"] = ogi
for i,v in enumerate(g.vertices()):
    ogi[v] = i

tmp = g.vertex(6) ## last vertex
g.vp.ogi[tmp] ## 6
g.vp["ogi"][tmp] ## 6
g.vp.ids[tmp] ## 'five'
g.vp.ids[6] ## 'five'
tmp2 = g.vertex('six') ## not valid function call; args must be coercable into an int
g.vp["ogi"][tmp2] ## 6
g.remove_vertex(6)
g.vertex_index[tmp] ## 6
tmp ## invalid vertex
tmp = g.vertex(4)
g.vp.ogi[tmp] ## 4
g.vp["ogi"][tmp] ## 4
g.vp.ids[tmp] ## 'three'
g.vp.ids[4] ## 'three'
g.vertex_index[tmp] ## 4
g.remove_vertex(3)
g.vertex_index[tmp] ## 4
tmp ## <Vertex object with index '4' at 0x7f8002b30640>


contract(g, e_property_map = econt_attr_map)

hl = g.new_vertex_property("bool")
g.vp["highlight"] = hl
to_hl = vertices(g, 4,5)
for v in g.vertices():
    hl[v] = v in to_hl

tmp = shrink(g, vertices(g, 4,5), v_property_map = {"leaf": lambda g,v:True}, e_property_map = ecoll_attr_map)
tmp.summarise_collapsed_vertices(g,"leaf")

## problem: we need to remove the collapsed vertices, but that results in loss of the vertex properties as well
## create a new object to store in collapsed_vertices where we copy over all the relevant properties?
g.remove_vertex(g.vertex(4))

gt.graphviz_draw(g, vcolor = g.vp["leaf"])


## read to graph_tool
import graph_tool.all as gt
g1 = gt.load_graph(graphml_1, fmt = "graphml")
g2 = gt.load_graph(graphml_2, fmt = "graphml")

g = gt.Graph([(1,2),(2,3),(3,1),(1,4)])
g_v2 = gt.Graph([(1,2),(2,3),(3,1),(1,4)])
g_v2.set_directed(False) ## undirected
gt.graphviz_draw(g) ## this one works. friendship ended with cairo, graphviz is my bestfriend now
gt.graph_draw(g) ## still doesn't work :(
# ## error message:
# TypeError: Couldn't find foreign struct converter for 'cairo.Context'
# Traceback (most recent call last):
#   File "/home/rachelle/miniconda3/envs/py3.13/lib/python3.13/site-packages/graph_tool/draw/gtk_draw.py", line 485, in layout_callback
#     self.regenerate_surface(reset=True, complete=True)
#     ~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/home/rachelle/miniconda3/envs/py3.13/lib/python3.13/site-packages/graph_tool/draw/gtk_draw.py", line 559, in regenerate_surface
#     self.base = w.create_similar_surface(cairo.CONTENT_COLOR_ALPHA,
#                 ~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
#                                          *geometry)
#                                          ^^^^^^^^^^
# TypeError: Couldn't find foreign struct converter for 'cairo.Surface'
# (<VertexPropertyMap object with value type 'vector<double>', for Graph 0x7f35e71347d0, at 0x7f35e71d9f40>, <VertexPropertyMap object with value type 'bool', for Graph 0x7f35e71347d0, at 0x7f35e71be580>)


## test how to check if is subgraph
g_tmp = gt.Graph([(1,2),(2,3)]) ## both edges are connected (by 2) and both exist in g
g_tmp2 = gt.Graph([(1,2),(2,4)]) ## (2,4) doesn't exist in g
g_tmp3 = gt.Graph([(2,3),(1,4)]) ## disjointed edges that both exist in g
g_tmp4 = gt.Graph([(1,2),(2,3),(1,4)])
g_tmp5 = gt.Graph([(2,1),(3,2)]) ## reverse direction
g_tmp6 = gt.Graph([(2,1),(3,2)])
g_tmp6.set_directed(False) ## undirected, but originally reverse direction
## Jaccard similarity = sum of equal non-zero entries in the adjacency matrices of both graphs
gt.similarity(g_tmp, g, asymmetric = False) ## 0.6666666666666666
## asymmstric = True --> asymmetric similarity of g1 to g2 will be computed (i.e. subgraph should be g1 to get 1.0)
gt.similarity(g_tmp, g, asymmetric = True) ## 1.0
gt.similarity(g_tmp3, g, asymmetric = True) ## 1.0; disjointed subgraph
gt.similarity(g_tmp5, g, asymmetric = True) ## 0.0; reverse direction
gt.similarity(g_tmp5, g_v2, asymmetric = True) ## 1.0; reverse direction; og tree undirected (but originally in opposite direction)
gt.similarity(g_tmp6, g, asymmetric = True) ## 0.5; subtree undirected (but edges originally reverse direction)
gt.similarity(g_tmp6, g_v2, asymmetric = True) ## 1.0; subtree and og tree undirected (but edges originally in opposite directions)

## adding properties to nodes/vertices
vaccid = g.new_vp("string")



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

def reduction_rule_1(f1, f2, disregard_root = True):
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


def Apx_MAF(f1, f2, disregard_root = True):
    """
    Input: two leaf-labelled forests F_1 and F_2 over the same label-set L
    Output: an L-partition P that is a solution for (F_1, F_2)
    1. apply Reduction Rules 1-2 and contractions on (F_1, F_2) until they
       are not applicable
    2. repeat until F_2 becomes a subgraph of F_1
    2.1 let X_1 be a BSS in F_1, and let X_2 be the leaf set in F_2 such that
        l(X_1) = l(X_2);
    2.2 switch
         for each Q, where Q = 1, 2.1, 2.2, 3, 4, listed in this order:
         for Case Q: apply Meta-Step Q;
    2.3 apply Reduction Rules 1-2 and contractions on (F_1, F_2) until they
        are not applicable
    3 return the L-partition P constructed from the label-partition for F_2
    """
    ## check if F_1 consists only of a single-vertex trees
    if all(t.is_single_vertex_tree() t for t in f1.trees):
        print("Unable to proceed: F_1 only consists of single-vertex trees")
        return None




###################
##  DEFINITIONS  ##
###################
