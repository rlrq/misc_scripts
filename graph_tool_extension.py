# from graph_tool.all import *

import graph_tool.all as gt

import math
import itertools
import functools
# import threading

# class Vertex(gt.Vertex):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

gt.Vertex.total_degree = lambda v:(v.in_degree() + v.out_degree())
# gt.Vertex.collapsed_vertices = []
# gt.Vertex.summarise_collapsed_vertices = lambda v,g,vprop:';'.join(str(g.vp[vprop][vcoll]) for vcoll in g.vp["collapsed_vertices"][v])

# class A():
#     def __init__(self, a):
#         self.a = a
#     def mult(self):
#         print("A.mult")
#         return self.a * 2

# class C():
#     def __init__(self, c):
#         self.c = c
#     def mult10(self):
#         print("C.mult10")
#         return self.c * 10
#     def mult(self):
#         print("C.mult")
#         return self.c * 5

# class B(A):
#     def __init__(self, b):
#         super().__init__(b)
#         self.c = C(b)
#     def mult(self):
#         print("B.mult")
#         return super().mult() * 2
#     def __getattr__(self, name):
#         return getattr(self.c, name)

# b = B(2)
# b.mult() ## B.mult; A.mult; 8
# A.mult(b) ## A.mult; 4
# b.mult10() ## C.mult10


def temporarily_undirected(func):
    def helper(*args, preserve_direction = False, **kwargs):
        ## get original directed-ness of all graphs passed as positional arguments
        if not preserve_direction:
            graphs = [a for a in args if isinstance(a, gt.Graph)]
            graphs_is_directed = {g: g.is_directed() for g in graphs}
            for g in graphs:
                g.set_directed(False)
        ## execute function
        try:
            output = func(*args, **kwargs)
            ## restore original directed-ness for all graphs passed as positional arguments
            if not preserve_direction:
                for g, directed in graphs_is_directed.items():
                    g.set_directed(directed)
            ## return function output
            return output
        ## restore original directedness before throwing error
        except Exception as e:
            ## restore original directed-ness for all graphs passed as positional arguments
            if not preserve_direction:
                for g, directed in graphs_is_directed.items():
                    g.set_directed(directed)
            raise e
    return helper

def invert_dict(d):
    invd = {}
    for k, v in d.items():
        invd[v] = invd.get(v, []) + [k]
    return invd

class FossilVertex():
    """Copies over all vertex properties so we retain a copy even after removing the original node"""
    def __init__(self, hv):
        self.hv = hv ## HashedVertex object (mostly obsolete once it's removed from the graph)
        self.hi = self.hv.hi ## the hashed id value
        self.hg = self.hv.hg ## HashedGraph object
        self.properties = {} ## {"<property name>": <property value>}
        self._copy_properties()
    def __getattr__(self, attr_name):
        """
        For any other attributes/methods/properties not explicitly defined in this class,
        check if attr_name is a vertex property (and return property value if so),
        else execute using gt.Vertex object.
        """
        ## for some reason, properties end up calling __getattr__ too even though this is only supposed
        ## to catch non-existent attributes
        if attr_name in dir(self): ## catch calls for class' properties
            return getattr(self, attr_name)
        elif attr_name in self.vp: ## redirect to vertex properties
            return self.vp[attr_name] ## retrieve the relevant value
        raise AttributeError(f"'{self.__class__.__name__}' has no attribute '{attr_name}'")
    @property
    def vp(self): return self.properties
    @property
    def vprop(self): return self.properties
    def _copy_properties(self):
        index = self.hv.index
        for vprop_name, vprop in self.hg.vp.items():
            self.properties[vprop_name] = vprop[index]
        return

@functools.total_ordering
class HashedVertex():
    ## This class was created so that the same object can always retrieve data of the same vertex regardless of re-indexing
    def __init__(self, hg, hashed_id):
        self.hg = hg ## HashedGraph object
        self.hi = hashed_id ## value that was hashed (ie stored in hg.vp.ids)
    def __eq__(self, other):
        return self.index == other.index
    def __lt__(self, other):
        return self.index < other.index
    def __hash__(self):
        return hash("HashedVertex_{id(self.hg)}_{self.hi}")
    def __repr__(self):
        return "<HashedVertex object with hashed id '%s' at 0x%x>" % (self.hi, id(self))
    def __getattr__(self, attr_name):
        """
        For any other attributes/methods/properties not explicitly defined in this class,
        check if attr_name is a vertex property (and return property value if so),
        else execute using gt.Vertex object.
        """
        ## for some reason, properties end up calling __getattr__ too even though this is only supposed
        ## to catch non-existent attributes
        if attr_name in dir(self): ## catch calls for class' properties
            return getattr(self, attr_name)
        elif attr_name in self.hg.vp: ## redirect to vertex properties in self.hg.vp
            return self.hg.vp[attr_name][self.index] ## retrieve the relevant value for self's vertex
        return getattr(self.v, attr_name) ## redirect to gt.Vertex object
    def __int__(self): ## return index if cast to integer
        return self.index
    @property
    def index(self):
        """Returns vertex index in graph"""
        # print(self.hi)
        return self.hg._hash_to_index(self.hi)
    @property
    def v(self):
        """Retrieves gt.Vertex object of self by calling Graph.vertex with correct index"""
        return self.hg.vertex(self.index)
    @property
    def graph(self): return self.hg
    @property
    def g(self): return self.hg
    def label(self, f):
        return f(self)

class CollapsedVertex(HashedVertex):
    def __init__(self, hg, hi, X, *args, **kwargs):
        super().__init__(hg, hi, *args, **kwargs)
        self._vertices = [FossilVertex(self.hg.hvertex(v)) for v in X]
    def __getattr__(self, attr_name):
        """
        For any other attributes/methods/properties not explicitly defined in this class,
        check if attr_name is a vertex property (and return property value if so),
        else execute using gt.Vertex object.
        """
        ## for some reason, properties end up calling __getattr__ too even though this is only supposed
        ## to catch non-existent attributes
        if attr_name in dir(self): ## catch calls for class' properties
            return getattr(self, attr_name)
        elif attr_name in self.hg.vp: ## redirect to vertex properties in self.hg.vp
            return self.hg.vp[attr_name][self.index] ## retrieve the relevant value for self's vertex
        return getattr(self.v, attr_name) ## redirect to gt.Vertex object
    def vertices(self):
        for v in self._vertices:
            if v.is_vcoll:
                for v2 in v.vcoll.vertices():
                    yield v2
            else:
                yield v
        return
    def vp_all(self, vprop_name):
        return [fv.vp.get(vprop_name, None) for fv in self.vertices()]
    def vp_flat(self, vprop_name):
        return self.vp_all(vprop_name)
    def label(self, label_f):
        return [label_f(v) for v in self.vertices()]

class HashedGraph(gt.Graph):
    def __init__(self, *args, hash_later = False, hash_leaves_only = False, relabel_internal = False, **kwargs):
        super().__init__(*args, hashed = (not hash_later), **kwargs)
        self._magnitude = math.ceil(math.log(len(self), 10))
        self.assign_leaf()
        self._internals_hashable = (not hash_leaves_only)
        if relabel_internal:
            ## relabels internal vertices uniquely so that they can be hashed
            self._relabel_internal_vertices()
        ## this variable decides whether to use original labels for hashing internal nodes (True = don't use)
        self.hash_leaves_only = hash_leaves_only
        self._hash = {} ## stores {<hashed id>: <index>}
        self._hash_hv = {} ## stores {<hashed id>: <HashedVertex>}
        self._num_vcoll = 0 ## stores number of CollapsedVertex objects (so we can give new ones a unique hash)
        if not hash_later:
            self._rehash()
            self._rehash_hv()
        self.threads = []
        self._add_vp_vcoll()
    def _relabel_internal_vertices(self, pattern = None):
        if pattern is None:
            pattern = "internal_{:0" + str(self._magnitude) + "d}"
        if "ids" not in self.vp:
            vprop_ids = self.new_vertex_property("string")
            self.vp["ids"] = vprop_ids
        else:
            vprop_ids = self.vp["ids"]
        internal_i = 0
        for v in self.vertices():
            if not self.vp.leaf[v]:
                vprop_ids[v] = pattern.format(internal_i)
                internal_i += 1
        self._internals_hashable = True
        return
    def set_vp_as_id(self, vp_name):
        if vp_name is None:
            return
        if "ids" not in self.vp:
            vprop_ids = self.new_vertex_property(vprop_source.value_type())
            self.vp["ids"] = vprop_ids
        else:
            vprop_ids = self.vp["ids"]
        vprop_source = self.vp[vp_name]
        for v in self.vertices():
            if self.hash_leaves_only and not self.vp.leaf[v]: continue
            vprop_ids[v] = vprop_source[v]
        self._rehash_all()
        return
    def _add_vp_vcoll(self): ## vcoll stores CollapsedVertex objects
        vcoll = self.new_vertex_property("object")
        self.vp["vcoll"] = vcoll
        is_vcoll = self.new_vertex_property("bool")
        self.vp["is_vcoll"] = is_vcoll
        return
    def _rehash_all(self):
        """Update both self._hash and self._hash_hv"""
        self._rehash()
        self._rehash_hv()
        return
    def _rehash(self):
        """Update self._hash"""
        for v in self.vertices():
            self._rehash_one_index(v)
        return
    def _rehash_one_index(self, v):
        """Update self._hash for a single vertex only. Takes gt.Vertex object."""
        # if self.hash_leaves_only and not self.vp.leaf[v]: return
        if (not self._internals_hashable) and (not self.vp.leaf[v]): retuern
        hash_id = self.vp.ids[v]
        if self._hash.get(hash_id, None) != int(v):
            self._hash[hash_id] = int(v)
        return
    def _rehash_hv(self):
        """Update self._hash_hv"""
        for v in self.vertices():
            self._rehash_one_hv(v)
        return
    def _rehash_one_hv(self, v):
        """Update self._hash_hv for a single vertex only. Takes gt.Vertex object."""
        # if self.hash_leaves_only and not self.vp.leaf[v]: return
        if (not self._internals_hashable) and (not self.vp.leaf[v]): return
        vertex_hash_id = self.vp.ids[v]
        if vertex_hash_id not in self._hash_hv:
            self._hash_hv[vertex_hash_id] = HashedVertex(self, vertex_hash_id)
        return
    def _check_hash_index(self, index, hi):
        """Checks if index points to vertex with a given hash id value (hi)"""
        return self.vp.ids[index] == hi
    def _hash_to_index(self, hi):
        """Retrieves vertex index by hashed id value"""
        index = self._hash.get(hi, None)
        ## check if index points to correct vertex, if not call _rehash to update hash table
        if not self._check_hash_index(index, hi):
            self._rehash()
            index = self._hash.get(hi, None) ## re-get index
        return index
    def _hash_to_hvertex(self, hi):
        """Retrieves HashedVertex object by hashed id value"""
        if hi not in self._hash_hv:
            self._rehash_hv()
        return self._hash_hv.get(hi, None)
    def assign_leaf(self):
        leaf = self.new_vertex_property("bool")
        for v in self.vertices():
            leaf[v] = v.out_degree() == 0 if self.is_directed() else v.total_degree() == 1
        self.vp["leaf"] = leaf
        return
    def hvertex(self, val):
        """
        Retrieve node by hashed id value or Vertex object. Returns a HashedVertex object.
        Returns val if val is already a HashedVertex object.
        Basically a retrieval/coercion function in one.
        """
        if isinstance(val, HashedVertex):
            return val
        elif isinstance(val, gt.Vertex):
            val = self.vp.ids[val]
        elif self._hash_to_index(val) is None:
            return None
        return self._hash_to_hvertex(val)
    def hvertices(self, leaves_only = False): ## analogous to gt.vertices
        for v in self.vertices():
            if leaves_only and not self.vp.leaf[v]: continue
            vertex_hash_id = self.vp.ids[v]
            ## hash HashedVertex if not hashed already
            if vertex_hash_id not in self._hash_hv:
                self._rehash_one_hv(v)
            yield self._hash_hv.get(vertex_hash_id, None)
        return
    def subset_hvertices(self, *vals):
        return [self.hvertex(val) for val in vals]
    def new_vcoll(self, X):
        vcoll_id = ("CollapsedVertex_{:0" + str(self._magnitude) + "d}").format(self._num_vcoll)
        vcoll = CollapsedVertex(self, vcoll_id, X)
        self._num_vcoll += 1
        return vcoll
    @temporarily_undirected
    def is_subgraph_v2(self, other, preserve_direction = False,
                       label_f1 = int, label_f2 = int, mask_internal = False):
        """
        Returns whether self is a subgraph of other by getting subgraph isomorphisms
        and then checking if the leaves match up.
        Compares vertices using vertex index by default.
        """
        ## get vertex labels (must be coercible into strings)
        svp1 = self.new_vertex_property("string")
        svp2 = other.new_vertex_property("string")
        for v in self.hvertices():
            svp1[v] = '' if (mask_internal and not self.vp.leaf[v]) else str(label_f1(v))
        for v in other.hvertices():
            svp2[v] = '' if (mask_internal and not other.vp.leaf[v]) else str(label_f2(v))
        ## get isomorphisms
        return len(gt.subgraph_isomorphism(self, other, vertex_label = (svp1, svp2))) > 0
    @temporarily_undirected
    def is_subgraph(self, other, preserve_direction = False,
                    label_f1 = int, label_f2 = int):
        """
        Returns whether self is a subgraph of other (i.e. asymmetric Jaccard similarity == 1).
        Compares edges using vertex index by default.
        """
        # g_is_directed = self.is_directed()
        # other_is_directed = other.is_directed()
        # if not preserve_direction:
        #     self.set_directed(False)
        #     other.set_directed(False)
        ## create scalar vertex property
        svp1 = self.new_vertex_property("int")
        svp2 = other.new_vertex_property("int")
        g_labels = [label_f1(hv) for hv in self.hvertices()]
        other_labels = [label_f2(hv) for hv in other.hvertices()]
        svp_map = {lab: i for i, lab in enumerate(g_labels)}
        for i in range(len(g_labels)):
            svp1[i] = i
        for i, lab in enumerate(other_labels):
            svp2[i] = svp_map[lab]
        ## check similarity
        jaccard_sim = gt.similarity(self, other, asymmetric = True, label1 = svp1, label2 = svp2)
        # self.set_directed(g_is_directed)
        # other.set_directed(other_is_directed)
        return jaccard_sim == 1
    def get_parent(self, X):
        if not self.is_sibling_set(X):
            # print("X is not a sibling set.")
            return None
        parents = list(set.intersection(*[set(v.all_neighbours()) for v in X]))
        if len(parents) != 1:
            print("Warning: X has more/less than one common connection.")
            return parents
        return parents[0]
    def all_bss(self):
        """Returns a list of all BSS [[HashedVertex objects in BSS 1], [HashedVertex objects in BSS 2]]"""
        ## get all leaves and their parents (leaves must only have one neighbour)
        single_edge_trees = []
        leaves = {}
        for hv in self.hvertices():
            ## skip vertex if not leaf
            if not hv.leaf: continue
            ## skip single vertex trees
            if hv.total_degree() == 0: continue
            parent = list(hv.all_neighbours())[0]
            if self.vp.leaf[parent]:
                single_edge_trees.append(tuple(sorted([int(hv), int(parent)])))
            leaves[int(hv)] = int(parent)
        ## group leaves into BSS
        bss = [tuple(sorted(v)) for v in invert_dict(leaves).values() if len(v) > 1]
        single_edge_trees = list(set(single_edge_trees))
        output = bss + single_edge_trees
        output = [self.subset_hvertices(*[self.vp.ids[v_i] for v_i in verts]) for verts in output]
        return output
    def remove_vertices(self, v_list):
        ## vertex removal updates index of every vertex after it
        ## so we shouldn't remove while iterating
        ## and for efficiency we should also remove from decreasing index order
        for v in reversed(sorted(v_list)):
            self.remove_vertex(v)
        return
    def orphan_vertex(self, v):
        """Remove all edges incident to a given vertex"""
        for e in v.all_edges():
            self.remove_edge(e)
        return
    @temporarily_undirected
    def in_same_subgraph(self, v_list):
        """Returns true if vertices in list are in same connected subgraph (i.e. component)"""
        vertex_indices = [self.vertex_index[v] for v in v_list]
        # g_is_directed = self.is_directed()
        # self.set_directed(False)
        components = gt.label_components(self)[0].a
        # self.set_directed(g_is_directed)
        return len(set(components[i] for i in vertex_indices)) == 1
    def vertices_of_single_edge_tree(self, v_list):
        """Returns true if vertices are from a single-edge tree. Assumes single-edge tree has exactly 1 edge and 2 vertices."""
        if not self.in_same_subgraph(v_list): return False
        if len(set(v_list)) != 2: return False
        neighbours = set(map(int, itertools.chain(*[v.all_neighbours() for v in v_list])))
        return neighbours == set(map(int, v_list))
    def is_sibling_set(self, X):
        """
        'A sibling set is a set of leaves that are all siblings.'
        'Two leaves in a forest are siblings if they either share the same parent, or are the two leaves of a single-edge tree.'
        """
        if self.vertices_of_single_edge_tree(X): return True ## return True if x's leaves are from the same single-edge tree
        if not self.in_same_subgraph(X): return False ## return False if vertices are from disjointed trees
        x = set(X)
        all_leaves = all(map(lambda v:(v.total_degree() == 1), x))
        if not all_leaves: return False
        neighbours = set(itertools.chain(*(v.all_neighbours() for v in x)))
        return len(neighbours) == 1 ## return whether siblings share the same parent
    def is_bss(self, X):
        """
        'A bottommost sibling set (abbr. BSS) is a maximal sibling set X such that either the degree of
          their parent is at most |X|+1, or X is the leaf set of a single-edge tree.'
        """
        ## return True if x's leaves are from the same single-edge tree
        if self.vertices_of_single_edge_tree(X): return True
        ## return False if vertices are from disjointed trees
        if not self.in_same_subgraph(X): return False
        X = set(X)
        all_leaves = all(map(lambda v:(v.total_degree() == 1), X))
        if not all_leaves: return False
        neighbours = set(itertools.chain(*[v.all_neighbours() for v in X]))
        if len(neighbours) != 1: return False ## return False if leaves don't all share the same parent
        return list(neighbours)[0].total_degree() <= (len(X)+1) ## return whether parent is at most |X|+1
    def shrink(self, X, v_property_map = {}, e_property_map = {}):
        """Replace vertices in X with a single vertex BUT ONLY IF X is a BSS."""
        if not self.is_bss(X): return
        parent = self.get_parent(X)
        ## create new vertex
        new_v = self.add_vertex()
        new_vcoll = self.new_vcoll(X)
        ## set vertex properties
        self.vp["vcoll"][new_v] = new_vcoll
        self.vp["ids"][new_v] = new_vcoll.hi
        self.vp["leaf"][new_v] = 1 ## new vertex is always a leaf because only BSSs are collapsed
        self.vp["is_vcoll"][new_v] = 1 ## new vertex is of CollapsedVertex class
        for name, f in v_property_map.items():
            self.vp[name][new_v] = f(self, new_v)
        ## connect new vertex to parent (if parent exists)
        if parent:
            new_e = self.add_edge(parent, new_v)
            ## set edge properties
            for name, f in e_property_map.items():
                self.ep[name][new_e] = f(self, new_e)
        ## remove vertices in X (connected edges will be automatically removed)
        self.remove_vertices(X)
        ## rehash all
        self._rehash_all()
        return self.hvertex(new_v) ## return the new vertex (as HashedVertex)
    def contract_vertex(self, v, e_property_map = {}):
        """
        Drop internal node with total degree 2, connect its neighbours with a new edge OR
        drop now-leaf nodes that weren't originally leaves.
        For undirected graphs only.
        Return number of contractions.
        """
        ## skip if leaf
        if self.vp["leaf"][v]:
            return 0
        ## skip if total degree != 2
        elif v.total_degree() > 2:
            return 0
        ## connect neighbouring vertices if total degree == 2
        elif v.total_degree() == 2:
            neighbours = list(v.all_neighbours())
            ## create new edge between v's neighbours
            new_e = self.add_edge(*neighbours)
            ## set edge properties
            for name, f in e_property_map.items():
                self.ep[name][new_e] = f(self, v)
        ## delete vertex (if total_degree <= 2)
        self.remove_vertex(v)
        return 1
    def contract(self, e_property_map = {}):
        """
        Drop internal nodes with total degree 2, connect its neighbours with a new edge.
        Drop now-leaf nodes that weren't originally leaves.
        Works properly for undirected graphs only.
        """
        num_contractions = 0
        ## vertex removal updates index of every vertex after it
        ## so we shouldn't remove while iterating
        ## and for efficiency we should also remove from decreasing index order
        v_list = reversed(sorted(self.vertices()))
        for v in v_list:
            num_contractions += self.contract_vertex(v, e_property_map = e_property_map)
        return num_contractions ## return the number of contractions so that we can determine when to stop contracting
    def draw(self, *args, device = "graphviz", background = False, penwidth = 2, vsize = 0.3,
             vcolour = None, ecolour = None, **kwargs):
        if device == "graphviz": dv = gt.graphviz_draw
        elif device == "cairo": dv = gt.cairo_draw
        elif device == "default": dv = gt.draw
        kwargs = {**kwargs, **{"penwidth": penwidth, "vsize": vsize}}
        if vcolour is not None:
            kwargs = {**{k: v for k, v in kwargs.items() if k != "vcolor"}, **{"vcolor": vcolour}}
        if ecolour is not None:
            kwargs = {**{k: v for k, v in kwargs.items() if k != "ecolor"}, **{"ecolor": vcolour}}
        if background: ## doesn't work. we still end up in graphviz's interface(?)
            t = threading.Thread(target = dv, name = device, args = (self,) + args, kwargs = kwargs)
            t.start()
            self.threads.append(t)
            return t
        else:
            dv(self, *args, **kwargs)
        return

class AxpMAFgraph(HashedGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return
    def meta_step_1(self, X):
        """Case 1. Leaves in X_2 are not in the same tree in F_2."""
        ## Let u_2 and v_2 be two leaves in X_2 that are in different trees in F_2,
        for i in range(len(X)-1):
            u = X[i]
            for j in range(i+1,len(X)):
                v = X[j]
                if not self.in_same_subgraph([u, v]):
                    ## then remove the edges incident to u_2 and v_2.
                    self.orphan_vertex(u)
                    self.orphan_vertex(v)
                    return 1
        return 0 ## operation not executed (yeah I know that 0 usually means no errors)
    def E_double_prime(self, X):
        """
        Sibling set X has parent p and the degree of p is at least |X|+2.
        Let E" be the set of edges that are incident to the parent p of X but not incident to the leaves in X.
        """
        if not self.is_sibling_set(X): return []
        parent = self.get_parent(X) ## 1 parent
        X_indices = set(int(v) for v in X)
        output = []
        for e in parent.all_edges():
            ## skip if edge is connected to any vertex in X
            if int(e.source()) in X_indices or int(e.target()) in X_indices: continue
            output.append(e)
        return output
    def meta_step_2(self, X):
        """Case 2. X_2 is a sibling set in F_2."""
        if (not self.is_bss(X)
            and self.is_sibling_set(X)):
            E_double_prime = self.E_double_prime(X)
        else:
            return 0 ## operation not executed
        ## get edges between X and parent
        X_indices = set(int(v) for v in X)
        parent = self.get_parent(X)
        E_X_parent = [e for e in parent.all_edges() if
                      (int(e.source()) in X_indices
                       or int(e.target()) in X_indices)]
        ## meta-step 2.1
        ## if |E"| > |X_2|, then pick any set E_1 of |X_2|-1 edges incident to the leaves in X_2,
        ## and pick any set E_2 of |X_2| edges in E", remove all edges in E_1 union E_2
        if len(E_double_prime) > len(X):
            E1 = E_X_parent[:-1] ## pick any set E_1 of |X_2|-1 edges incident to leaves in X
            E2 = E_double_prime[:len(X)] ## pick any set E_2 of |X_2| edges in E"
            for e in E1 + E2:
                self.remove_edge(e)
        ## meta-step 2.2
        ## if |E"| <= |X_2|, then pick any set E'_1 of |E"|-1 edges incident to the leaves in X_2,
        ## and remove all edges in E'_1 union E".
        else:
            E_prime_1 = E_X_parent[:(len(E_double_prime)-1)] ## pick any set E'_1 of |E"|-1 edges incident to the leaves in X_2
            for e in E_prime_1 + E_double_prime:
                self.remove_edge(e)
        return 1 ## operation executed successfully
    @temporarily_undirected
    def meta_step_3(self, X, preserve_direction = False):
        """Case 3. X_2 contains a sibling set Y_2 with |Y_2|>=2. Assumption: X_2 is not a sibling set."""
        ## Assumption: X is not a sibling set. assume also: graph is acyclic
        ## Find a sibling set
        ## 'Let u_2, v_2 be elements of Y_2'
        if len(X) < 3: return 0 ## operation not executed. Too few vertices
        Y = None
        break_loop = False
        for i in range(len(X)-1):
            u = X[i]
            for j in range(i+1,len(X)):
                v = X[j]
                if self.is_sibling_set([u, v]):
                    Y = [u, v]
                    break_loop = True
                    break
            if break_loop: break
        if Y is None: return 0 ## operation not executed. Premise that X_2 contains a sibling set Y_2 is false
        ## 'Let w_2 be a leaf in X_2 that is not a sibling of u_2 and v_2.'
        Y_parent = self.get_parent(Y) ## 1 parent
        for w in X:
            w_neighbour_indices = [int(n) for n in w.all_neighbours()]
            if Y_parent not in w_neighbour_indices: break ## get first vertex without u+v's parent as a neighbour
        ## 'Let e be an edge incident to the parent of w_2 but not on the path between u_w and w_2.'
        # g_is_directed = self.is_directed()
        # self.set_directed(False) ## set undirected to simplify steps
        path_vertices = gt.shortest_path(self, u, w)[0]
        path_vertices_index = [int(n) for n in path_vertices]
        w_parent = path_vertices[-2] ## -1 is w
        for n in w_parent.all_neighbours():
            ## check if n is along path between u & w
            if int(n) not in path_vertices_index:
                ## if not along path, get edge between n & w_parent
                e = self.edge(n, w_parent)
                break
        ## 'Then remove the edge e, the edge e_u incident to u, and the edge e_w incident to w'
        self.remove_edge(e)
        self.remove_edge(self.edge(u, Y_parent))
        self.remove_edge(self.edge(w, w_parent))
        # ## set direction back to original (even though this operation only really works on undirected graphs anyway...)
        # self.set_directed(g_is_directed)
        return 1 ## operation executed successfully
    @temporarily_undirected
    def meta_step_4(self, X, preserve_direction = False):
        """Case 4. All leaves in X_2 are in the same tree in F_2 but no two are siblings"""
        ## Assumption: there are no sibling sets.
        ## Let u_2, v_2 be a elements of X_2 such that the parent p_2 of u_2 has degree 2 in F_2[l(X_2)].
        ## (l(X_2) is the set of labels for leaves X_2)
        ## (F_2[l(X_2)] is presumably the tree prune of all leaves except X_2; subforest of F_2 induced by l(X_2))
        for i in range(len(X)-1):
            u = X[i]
            for j in range(i+1,len(X)):
                v = X[j]
                if self.is_sibling_set([u, v]):
                    return 0 ## sibling sets violate the premise
        ## Let e be an edge in F_2 that is incident to p_2 but not on the path P_{uv} between u_2 and v_2.
        # g_is_directed = self.is_directed()
        # self.set_directed(False) ## set undirected to simplify steps
        path_vertices = gt.shortest_path(self, u, v)[0]
        path_vertices_index = [int(n) for n in path_vertices]
        p = path_vertices[1] ## 0 is u, -1 is v
        e = None
        for n in p.all_neighbours():
            ## check if n is along path between u & w
            if int(n) not in path_vertices_index:
                ## if not along path, get edge between n & p
                e = self.edge(n, p)
                break
        ## Then cut the edge e, the edge e_u incident to u_2, and the edge e_v incident to v_2.
        if e is None: return 0
        # try:
        #     self.remove_edge(e)
        # except Exception as er:
        #     print(path_vertices)
        #     print(p.all_neighbours())
        #     print([int(n) for n in path_vertices])
        #     print(p)
        #     print(n)
        #     raise er
        self.remove_edge(self.edge(u, p))
        self.remove_edge(self.edge(v, path_vertices[-2]))
        # ## set direction back to original (even though this operation only really works on undirected graphs anyway...)
        # self.set_directed(g_is_directed)
        return 1 ## operation executed successfully
    def meta_step_switch(self, X, preserve_direction = False):
        """Switch for meta-steps"""
        executed = self.meta_step_1(X) ## if leaves are not in same tree
        if not executed: executed = self.meta_step_2(X) ## if X_2 is a sibling set
        if not executed: executed = self.meta_step_3(X) ## if X_2 contains a sibling set
        if not executed: executed = self.meta_step_4(X) ## if all leaves are in same tree but no 2 are siblings
        ## return whether any meta-steps were executed (0/1)
        return executed

def label_as_list(hv, label_f):
    """
    Takes a HashedVertex (hv) and a function that generates a label (label_f).
    Returns label(s) as a list.
    """
    if hv.is_vcoll:
        return hv.vcoll.label(label_f)
    return [hv.label(label_f)]

def flatten_hvertices(X, label_f):
    """
    Takes a list of HashedVertices (X) and a function that generates a label (label_f).
    Expands labels of collapsed vertices to the same level as regular vertices.
    """
    output = []
    for hv in X:
        output.extend(label_as_list(hv, label_f))
    return output

def label_is_subset(hv, label_f, labels):
    """
    Takes a HashedVertex, a function that generates a label (label_f), and an iterable of labels (labels).
    Checks if the vertex's label(s) is a subset of the list of labels.
    """
    hv_labels = set(label_as_list(hv, label_f))
    return set(hv_labels).issubset(set(labels))

## Forests F1, F2 must be HashedGraph objects
def reduction_rule_1(F1, F2, label_f1 = int, label_f2 = int):
    """
      If a label l is in a single-vertex tree in one of the forests F_1 and F_2,
    then remove the edge (if any) incident to the label l in the other forest.
    """
    ## get labels from single-vertex trees
    def single_vertex_labels(F, label_f):
        labels = []
        for hv in F.hvertices():
            if hv.total_degree() == 0:
                labels.extend(label_as_list(hv, label_f))
        return labels
    ## delete all edges adjacent to vertices whose label value correspond to single vertex trees in the other forest
    def prune(F, label_f, labels):
        labels = set(labels)
        edges_removed = 0
        for hv in F.hvertices():
            if label_is_subset(hv, label_f, labels):
                edges = reversed(sorted(hv.all_edges()))
                for e in edges:
                    F.remove_edge(e)
                    edges_removed += 1
        return edges_removed ## return this so we can determine when to stop
    ## reduce repeatedly until no more reduction is possible
    F1_single_vertex_labels = single_vertex_labels(F1, label_f1)
    F2_single_vertex_labels = single_vertex_labels(F2, label_f2)
    p1 = prune(F1, label_f1, F2_single_vertex_labels)
    p2 = prune(F2, label_f2, F1_single_vertex_labels)
    return p1 + p2 ## return this so we know if any reduction happened

## apply reduction rule 1 until it no longer applies
def reduction_rule_1_exhaustive(F1, F2, **kwargs):
    edges_removed = 0
    while True:
        newly_removed = reduction_rule_1(F1, F2, **kwargs)
        edges_removed += newly_removed
        if newly_removed == 0: break
    return edges_removed

def get_X2(F1, F2, X1, label_f1 = int, label_f2 = int):
    labels_X1 = set(flatten_hvertices(X1, label_f1))
    X2 = [hv for hv in F2.hvertices() if label_is_subset(hv, label_f2, labels_X1)]
    return X2

## Forests F1, F2 must be HashedGraph objects
def reduction_rule_2(F1, F2, X1, label_f1 = int, label_f2 = int,
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
      X1 is a reusable iterable of HashedVertex objects.
    """
    if not F1.is_bss(X1):
        print("X1 is not a BSS.")
        print([F1.vp.ids[v] for v in X1])
        return 0 ## operation not executed
    ## get equivalent X2 vertices
    X2 = get_X2(F1, F2, X1, label_f1 = label_f1, label_f2 = label_f2)
    if not F2.is_bss(X2):
        # print("X2 is not a BSS.") ## don't print this. There may be too many error messages.
        return 0 ## operation not executed
    ## shrink (HashedGraph.shrink returns a HashedVertex object)
    vcoll1 = F1.shrink(X1, v_property_map = vcoll_property_map, e_property_map = ecoll_property_map)
    vcoll2 = F2.shrink(X2, v_property_map = vcoll_property_map, e_property_map = ecoll_property_map)
    ## contract (if necessary)
    v1 = None if len(v:=list(vcoll1.all_neighbours())) == 0 else v[0]
    v2 = None if len(v:=list(vcoll2.all_neighbours())) == 0 else v[0]
    if v1: _ = F1.contract_vertex(v1, e_property_map = econt_property_map)
    if v2: _ = F2.contract_vertex(v2, e_property_map = econt_property_map)
    return 1

## apply reduction rule 2 on all BSS in F1
def reduction_rule_2_all_bss(F1, F2, **kwargs):
    bss_shrunk = 0
    bss_F1 = F1.all_bss()
    if len(bss_F1) == 0:
        # print("No BSS in F1.")
        return 0
    for X1 in bss_F1:
        bss_shrunk += reduction_rule_2(F1, F2, X1, **kwargs)
    return bss_shrunk

## apply reduction rule 2 on all BSS in F1 repeatedly until it no longer applies
def reduction_rule_2_exhaustive(F1, F2, **kwargs):
    bss_shrunk = 0
    while True:
        newly_shrunk = reduction_rule_2_all_bss(F1, F2, **kwargs)
        bss_shrunk += newly_shrunk
        if newly_shrunk == 0: break
    return bss_shrunk

def apply_reduction_rules_exhaustive(F1, F2, label_f1 = int, label_f2 = int,
                                     vcoll_property_map = {}, ecoll_property_map = {},
                                     econt_property_map = {}):
    cycles = 0
    total_removed = 0
    total_shrunk = 0
    while True:
        edges_removed = reduction_rule_1_exhaustive(F1, F2, label_f1 = label_f1, label_f2 = label_f2)
        bss_shrunk = reduction_rule_2_exhaustive(F1, F2, label_f1 = label_f1, label_f2 = label_f2,
                                                 vcoll_property_map = vcoll_property_map,
                                                 ecoll_property_map = ecoll_property_map,
                                                 econt_property_map = econt_property_map)
        cycles += 1
        total_removed += edges_removed
        total_shrunk += bss_shrunk
        if edges_removed == bss_shrunk == 0: break
    return cycles, total_removed, total_shrunk

def apply_meta_steps_on_bss(F1, F2, label_f1 = int, label_f2 = int,
                            econt_property_map = {}, **rule_2_kwargs):
    bss_iter = 0
    executed = 0
    total_removed = 0
    total_shrunk = 0
    for X1 in F1.all_bss():
        bss_iter += 1
        X2 = get_X2(F1, F2, X1, label_f1 = label_f1, label_f2 = label_f2)
        executed = F2.meta_step_switch(X2)
        if executed: break
    ## apply reduction rules 1-2 and contractions on (F1, F2) until they are no longer applicable
    F1.contract(e_property_map = econt_property_map)
    F2.contract(e_property_map = econt_property_map)
    cycles, removed, shrunk = apply_reduction_rules_exhaustive(F1, F2, label_f1 = label_f1, label_f2 = label_f2,
                                                               econt_property_map = econt_property_map, **rule_2_kwargs)
    total_removed += removed
    total_shrunk += shrunk
    return bss_iter, executed, total_removed, total_shrunk

def Apx_MAF(F1, F2, label_f1 = int, label_f2 = int, internal_unlabelled = True,
            econt_property_map = {}, **rule_2_kwargs):
    ## check if F1 and F2 have the same label-set
    if set(label_f1(hv) for hv in F1.hvertices() if F1.vp.leaf[hv]) != set(label_f2(hv) for hv in F2.hvertices() if F2.vp.leaf[hv]):
        print("F1 and F2 do NOT have the same label-set.")
        return
    label_f_kwargs = {"label_f1": label_f1, "label_f2": label_f2}
    ## apply reduction rules 1-2 and contractions on (F1, F2) until they are no longer applicable
    F1.contract(e_property_map = econt_property_map)
    F2.contract(e_property_map = econt_property_map)
    _, total_removed, total_shrunk = apply_reduction_rules_exhaustive(F1, F2, **label_f_kwargs, **rule_2_kwargs)
    total_operations = 0
    ## execute meta-steps and reduction rules until F2 is a subgraph of F1
    while not F2.is_subgraph_v2(F1, mask_internal = internal_unlabelled, **label_f_kwargs):
        bss_iter, executed, removed, shrunk = apply_meta_steps_on_bss(F1, F2, **label_f_kwargs)
        total_operations += executed
        total_removed += removed
        total_shrunk += shrunk
        if not bss_iter or not executed:
            break
    ## get L-partitions (label partitions)
    F2_is_directed = F2.is_directed()
    F2.set_directed(False)
    components = gt.label_components(F2)[0].a
    F2.set_directed(F2_is_directed)
    L_partition = {}
    for v_i, component in enumerate(components):
        L_partition[component] = L_partition.get(component, []) + label_as_list(F2.hvertex(F2.vertex(v_i)), label_f2)
    L_partition = list(L_partition.values())
    return L_partition, total_operations, total_removed, total_shrunk

# ## some tests
# g = HashedGraph([("one","two"),("one","three"),("one","nine"),("zero","one"),("zero","four"),("four","five"),("four","six"),("nine","seven"),("seven","eight"),("seven","ten")])
# g.assign_leaf()
# g.vp.leaf[g.hvertex("two").index]
# g.vp.leaf[g.hvertex("two")] ## we can do this because HashedVertex.__int__ returns the index
# # x = g.draw()

# ogi = g.new_vertex_property("int")
# g.vp["ogi"] = ogi
# for i,v in enumerate(g.vertices()):
#     ogi[v] = i

# # vcoll = g.new_vertex_property("object")
# # g.vp["vcoll"] = vcoll

# x = g.hvertex("six")
# g.vp.ogi[x] ## 7
# x.ogi ## 7 (redirected by __getattr__ to g.vp.ogi[x])
# x.total_degree() ## correctly calls gt.Vertex method (well, a custom method that we added to it)
# g.is_bss([g.hvertex("eight"), g.hvertex("seven")])
# tmp = g.shrink([g.hvertex("eight"), g.hvertex("seven")]) ## should return None, because "eight" and "seven" isn't a bss
# tmp = g.shrink([g.hvertex("eight"), g.hvertex("ten")]) ## should return a HashedVertex object
# tmp.hi ## 'CollapsedVertex_00'
# tmp.vcoll.vp_all("leaf") ## [1,1] returns a list of the values of a given property of the vertices it collapsed
# # g.draw(vcolor = g.vp.leaf)

# g.contract() ## 3
# # g.draw(vcolor = g.vp.is_vcoll)

# g.is_bss([g.hvertex("two"), g.hvertex("three"), g.hvertex("CollapsedVertex_00")])
# tmp2 = g.shrink([g.hvertex("two"), g.hvertex("three"), g.hvertex("CollapsedVertex_00")])
# tmp2.vcoll.vp_all("ids") ## ['two', 'three', 'eight', 'ten'] successfully expands all collapsed vertices (even 2nd level ones)
# # g.draw(vcolor = g.vp.is_vcoll)
# # g.draw(vcolor = g.vp.leaf)

# g1 = HashedGraph([("one","two"),("one","three"),("one","nine"),("zero","one"),("zero","four"),("four","five"),("four","six"),("nine","seven"),("seven","eight"),("seven","ten")])
# g2 = HashedGraph([("zero","one"),("zero","four"),("four","five"),("four","six"),("nine","seven"),("seven","eight"),("seven","ten"),("one","two"),("one","three"),("one","nine")]) ## different edge/vertex order
# g3 = HashedGraph([("zero","one"),("zero","four"),("four","five"),("four","six"),("nine","seven"),("seven","eight"),("seven","ten"),("one","two"),("one","three"),("nine","one")]) ## 1 different edge direction from g2
# g2.is_subgraph(g1) ## False; default comparison is between vertex indices
# g1.is_subgraph(g2, label_f1 = lambda hv:hv.hi, label_f2 = lambda hv:hv.hi) ## True
# g2.is_subgraph(g1, label_f1 = lambda hv:hv.hi, label_f2 = lambda hv:hv.hi) ## True
# g3.is_subgraph(g2) ## True (default behaviour compares undirected versions of the graphs)
# g3.is_subgraph(g2, preserve_direction = True) ## False
# g2.remove_edge(g2.edge(g2.hvertex("one"), g2.hvertex("two"))) ## now edges g2 < g1
# g2.is_subgraph(g1, label_f1 = lambda hv:hv.hi, label_f2 = lambda hv:hv.hi) ## True
# g1.is_subgraph(g2, label_f1 = lambda hv:hv.hi, label_f2 = lambda hv:hv.hi) ## False

# label_no_internals = lambda hv:('' if not hv.leaf else hv.ids)
# ## test how providing some labels work
# ## I like it. It matches the vertices to the labels and if multiple have the same label then only their relative toplogy to the vertices with unique labels matter.
# l2 = lambda hv:('' if not hv.leaf else (hv.ids if hv.ids in ["two", "three"] else '')) ## 4 maps
# l2 = lambda hv:('' if not hv.leaf else (hv.ids if hv.ids in ["two", "three", "five"] else '')) ## 2 maps
# svp1 = g1.new_vertex_property("string")
# for v in g1.hvertices():
#     svp1[v] = str(l2(v))

# svp2 = g2.new_vertex_property("string")
# for v in g2.hvertices():
#     svp2[v] = str(l2(v))

# svp3 = g3.new_vertex_property("string")
# for v in g3.hvertices():
#     svp3[v] = str(l2(v))

# gt.subgraph_isomorphism(g1, g2, vertex_label = (svp1, svp2))
# gt.subgraph_isomorphism(g1, g3, vertex_label = (svp1, svp3)) ## no map
# g1.set_directed(False); g3.set_directed(False)
# gt.subgraph_isomorphism(g1, g3, vertex_label = (svp1, svp3)) ## 2 maps


## if vertex is a collapsed vertex, use vertex's hashed ID instead of applying label_f
def use_vcoll_id(label_f):
    def helper(hv):
        if hv.is_vcoll:
            return hv.hi
        return label_f(hv)
    return helper

def make_in_out_scalar(g, X):
    vp = g.new_vertex_property("bool")
    X_indices = set(map(int, X))
    for v in g.vertices():
        vp[int(v)] = int(v) in X_indices
    return vp

def make_groups_scalar(g, *groups):
    groups_inv = {int(v): i+1 for i, verts in enumerate(groups) for v in verts}
    vp = g.new_vertex_property("int")
    for v in g.vertices():
        vp[int(v)] = groups_inv.get(int(v), 0)
    return vp

# ## test meta steps
# zero = "zero"; one = "one"; two = "two"; three = "three"; four = "four"; five = "five"; six = "six"; seven = "seven"; eight = "eight"; nine = "nine"; ten = "ten"; eleven = "eleven"; twelve = "twelve"; thirteen = "thirteen"; fourteen = "fourteen"
# g = HashedGraph([(zero,one),(zero,two),(zero,three),(zero,four),(zero,five),(zero,six),(six,seven),(six,eight),(six,nine),(six,ten),(ten,eleven),(ten,twelve),(ten,thirteen)])
# X2 = g.subset_hvertices(one, thirteen)
# vp_groups = make_groups_scalar(g, X2, g.subset_hvertices(two, three, four, five), g.subset_hvertices(seven, eight, nine), g.subset_hvertices(eleven, twelve))
# g.meta_step_2(X2)
# g.meta_step_3(X2)
# g.meta_step_4(X2)
# # g.draw(vcolor = vp_groups)

# g2 = HashedGraph([(zero,one),(zero,two),(zero,three),(zero,four),(zero,five),(zero,six),(six,seven),(six,eight),(six,nine),(six,ten),(ten,eleven),(ten,twelve),(ten,thirteen)])
# reduction_rule_1(g, g2, label_f1 = lambda hv:hv.hi, label_f2 = lambda hv:hv.hi)

# ## verify that we're shrinking the correct bss
# g.draw(vcolor = make_groups_scalar(g, g.subset_hvertices(eleven, twelve)))
# ## shrink
# tmp = g.shrink(g.subset_hvertices(eleven, twelve))
# # g.draw(vcolor = g.vp.leaf)
# g.contract()
# # g.draw(vcolor = g.vp.leaf)
# g.draw(vcolor = make_groups_scalar(g, *map(lambda x:[x], g.subset_hvertices(zero, one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve))))




#### for development, use these trees, which have some multifurcating subtrees
import sys
sys.path.append("/mnt/chaelab/rachelle/src")
import newick_extension as ne
import newick
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

## convert to edge list (text file, format based on https://docs.oracle.com/en/database/oracle/property-graph/23.1/spgdg/edge-list-edge_list.html)
def parse_if_non_numeric(v, f):
    try:
        v = float(v)
        return ''
    except ValueError:
        return f(v)
    except TypeError:
        return ''

def edgelist_to_hashed_graph(fin, hash_id, graph_class = HashedGraph,
                             hash_leaves_only = True, relabel_internal = True, **kwargs):
    v_properties, e_properties, v_dict, e_list = ne.parse_edgelist(fin, **kwargs)
    eprops = [(eprop, ne.reformat_type(ne.infer_type(e[i+2] for e in e_list).__name__))
              for i, eprop in enumerate(e_properties)]
    hash_id_index = v_properties.index(hash_id)
    g = graph_class(e_list, eprops = eprops, hash_later = True,
                    hash_leaves_only = hash_leaves_only, relabel_internal = relabel_internal)
    ## assign vertex properties
    ogi = g.new_vertex_property("int") ## original index
    g.vp["ogi"] = ogi
    for v in g.vertices():
        ogi[v] = int(v)
    ## properties from edgelist data
    for vprop_i, v_prop in enumerate(v_properties):
        ## create property
        vprop = g.new_vertex_property(
            ne.reformat_type(ne.infer_type(v[vprop_i] for v in v_dict.values()).__name__))
        g.vp[v_prop] = vprop
        ## iterate through vertices and assign property values
        for v_i, v_dat in v_dict.items():
            vert = g.vertex(v_i) ## get vertex
            vprop[vert] = v_dat[vprop_i] ## assign
    g.set_vp_as_id(hash_id)
    return g


v_attr_map = {
    "seqlen": lambda n:seqlens.get(n.name, ''),
    "accid": lambda n:parse_if_non_numeric(n.name, get_accid),
    "gid": lambda n:parse_if_non_numeric(n.name, get_gid),
    "leaf": lambda n:n.is_leaf# ,
    # "collapsed_vertices": lambda n:[]
}
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
    "leaf": lambda g,v:True# ,
    # "collapsed_vertices": lambda g,n:[]
}
ecoll_attr_map = {
    "length": lambda g,e:-1
}
econt_attr_map = {
    "length": lambda g,v:sum(g.ep["length"][e] for e in v.all_edges())
}

## read & process trees
t1 = newick.read(nwk1)[0]
t1.set_leaf_attributes(attr_map)
t1.prune_leaves_by_group("accid", lambda leaves:max(leaves, key = lambda leaf:leaf.seqlen), in_place = True)
t1.match = lambda n:n.accid
t2 = newick.read(nwk2)[0]
t2.set_leaf_attributes(attr_map)
t2.prune_leaves_by_group("accid", lambda leaves:max(leaves, key = lambda leaf:leaf.seqlen), in_place = True)
ne.prune_unique_leaves(t1, t2, match_function = lambda n:n.accid, in_place = True)
t2.match = lambda n:n.accid
## let's try something, where we use the same tree for t1 and t2 but swap two leaf labels in t2
t3 = newick.read(nwk1)[0]
t3.set_leaf_attributes(attr_map)
t3.prune_leaves_by_group("accid", lambda leaves:max(leaves, key = lambda leaf:leaf.seqlen), in_place = True)
ne.prune_unique_leaves(t1, t3, match_function = lambda n:n.accid, in_place = True)
t3.match = lambda n:n.accid
t3_leaves = t3.get_leaves()
t1_leaves = t1.get_leaves()
t3_leaves[0].name = t1_leaves[20].name ## 'Reference|ANGE-B-10.g40|CDS|ANGE-B-10.g40.t1'
t3_leaves[20].name = t1_leaves[0].name ## 'Reference|SALE-A-10.g39|CDS|SALE-A-10.g39.t1'

## convert to edgelist
tmp1 = "/mnt/chaelab/rachelle/tmp1.txt"
tmp2 = "/mnt/chaelab/rachelle/tmp2.nwk"
tmp3 = "/mnt/chaelab/rachelle/tmp3.nwk"
edgelist_1 = "/mnt/chaelab/rachelle/tmp/t1.edgelist"
edgelist_2 = "/mnt/chaelab/rachelle/tmp/t2.edgelist"
edgelist_3 = "/mnt/chaelab/rachelle/tmp/t3.edgelist"
newick.write(t1, tmp1)
newick.write(t2, tmp2)
newick.write(t3, tmp3)
ne.newick_to_edgelist(tmp1, edgelist_1, **edgelist_kwargs)
ne.newick_to_edgelist(tmp2, edgelist_2, **edgelist_kwargs)
ne.newick_to_edgelist(tmp3, edgelist_3, **edgelist_kwargs)

## read edgelist into graph-tool graph
g1 = edgelist_to_hashed_graph(edgelist_1, "label", graph_class = AxpMAFgraph,
                              v_property_coerce = v_property_coerce,
                              e_property_coerce = e_property_coerce)
g2 = edgelist_to_hashed_graph(edgelist_2, "label", graph_class = AxpMAFgraph,
                              v_property_coerce = v_property_coerce,
                              e_property_coerce = e_property_coerce)
g3 = edgelist_to_hashed_graph(edgelist_3, "label", graph_class = AxpMAFgraph,
                              v_property_coerce = v_property_coerce,
                              e_property_coerce = e_property_coerce)
g0 = edgelist_to_hashed_graph(edgelist_1, "label", graph_class = AxpMAFgraph,
                              v_property_coerce = v_property_coerce,
                              e_property_coerce = e_property_coerce)

## undirected
g1.set_directed(False)
g2.set_directed(False)
g3.set_directed(False)

# ## draw
# g1.draw()
# g2.draw()

## the pair of trees with one pair of leaves swapped
am_1swap = Apx_MAF(g1, g3, label_f1 = use_vcoll_id(lambda hv:hv.accid), label_f2 = use_vcoll_id(lambda hv:hv.accid))
am_1swap[1:] ## (2, 5, 68) ## L-partitions, # meta step operations, # single-vertex trees removed, # bss shrunk
g1.draw()
g3.draw()
g0.draw(vcolour = make_groups_scalar(g0, g0.subset_hvertices('Reference|ANGE-B-10.g40|CDS|ANGE-B-10.g40.t1', 'Reference|SALE-A-10.g39|CDS|SALE-A-10.g39.t1')))
tmp2 = [(g0.subset_hvertices(v.hi) if not v.is_vcoll else g0.subset_hvertices(*v.vcoll.vp_all("ids"))) for v in g1.hvertices()]
g0.draw(vcolour = make_groups_scalar(g0, *tmp2))


## the pair of trees from different ADR1 genes
g4 = edgelist_to_hashed_graph(edgelist_1, "label", graph_class = AxpMAFgraph,
                              v_property_coerce = v_property_coerce,
                              e_property_coerce = e_property_coerce)
g4.set_directed(False)
am_adr1s = Apx_MAF(g4, g2, label_f1 = use_vcoll_id(lambda hv:hv.accid), label_f2 = use_vcoll_id(lambda hv:hv.accid))
am_adr1s[1:] ## (54, 124, 1) ## L-partitions, # meta step operations, # single-vertex trees removed, # bss shrunk
g4.draw()
g2.draw()
tmp2 = [(g0.subset_hvertices(v.hi) if not v.is_vcoll else g0.subset_hvertices(*v.vcoll.vp_all("ids"))) for v in g4.hvertices()]
g0.draw(vcolour = make_groups_scalar(g0, *tmp2))


# g2_comp = gt.label_components(g2)[0].a
# g2_comp_d = {g2.vp.accid[i]: int(v) for i, v in enumerate(g2_comp) if g2.vp.leaf[i]}
# g2_comp_d_rev = invert_dict(g2_comp_d)
# g1_accid_map = {g1.vp.accid[v]: g1.vp.ids[v] for v in g1.vertices() if g1.vp.leaf[v]}
# g2_comp_g1 = {comp_i: [g1_accid_map[accid] for accid in accids] for comp_i, accids in g2_comp_d_rev.items()}

# g1.draw(vcolor = make_groups_scalar(g1, g1.subset_hvertices(*g2_comp_g1[0]), g1.subset_hvertices(*g2_comp_g1[1]), g1.subset_hvertices(*g2_comp_g1[3])))
# g1.draw(vcolor = g1.vp.leaf)
# g2.draw(vcolor = g2.vp.leaf)
# g1.draw(vcolor = g1.vp.is_vcoll)
# g2.draw(vcolor = g2.vp.is_vcoll)

# label_f_kwargs = {"label_f1": use_vcoll_id(lambda hv:hv.accid), "label_f2": use_vcoll_id(lambda hv:hv.accid)}
# apply_reduction_rules_exhaustive(g1, g2, **label_f_kwargs)



class ParaMAFgraph(HashedGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return
    def cut(self, hv):
        pass
    def branching_rule_1(self, X):
        """Case 1. Leaves in X_2 are not in the same tree in F_2."""
        ## Let u_2 and v_2 be two leaves in X_2 that are in different trees in F_2,
        for i in range(len(X)-1):
            u = X[i]
            for j in range(i+1,len(X)):
                v = X[j]
                if not self.in_same_subgraph([u, v]):
                    ## then decrease k by 1, and branch two ways: [W1] cut u_2 in F_2; and [W2] cut v_2 in F_2.
                    pass
                    return 1
        return 0 ## operation not executed (yeah I know that 0 usually means no errors)
    def branching_rule_2(self, X):
        """Case 2. X_2 is a sibling set in F_2."""
        if (not self.is_bss(X)
            and self.is_sibling_set(X)):
            E_double_prime = self.E_double_prime(X)
        else:
            return 0 ## operation not executed
        ## fix a leaf u_2 in X_2, and let E_2 = {e_1,...,e_h}, h>=2, be the set of edges that are incident to the
        ## parent p_2 of X_2 but not to any leaf in X_2. Then branch into h+1 ways:
        ## [W(0)] remove all edges incident to X_2 except the one incident to u_2 and decrease k by |X_2|-1;
        ## and, for each 1<=i<=h, [W(i)] remove all edges in E_2 except e_i and decrease k by h-1
        #### TODO
        return 1
    def branching_rule_3(self, X):
        """Case 3. X_2 contains a sibling set Y_2 with |Y_2|>=2. Assumption: X_2 is not a sibling set."""
        ## Assumption: X is not a sibling set. assume also: graph is acyclic
        ## Find a sibling set
        ## 'Let u_2, v_2 be elements of Y_2'
        if len(X) < 3: return 0 ## operation not executed. Too few vertices
        Y = None
        break_loop = False
        for i in range(len(X)-1):
            u = X[i]
            for j in range(i+1,len(X)):
                v = X[j]
                if self.is_sibling_set([u, v]):
                    Y = [u, v]
                    break_loop = True
                    break
            if break_loop: break
        if Y is None: return 0 ## operation not executed. Premise that X_2 contains a sibling set Y_2 is false
        #### TODO: validate that the bit of code above is in line with the algorithm
        ## branch into three ways:
        ## [W1] cut all leaves in Y_2 except u_2 in F_2 and decrease k by |Y_2|-1;
        ## [W2] cut z_2 in F_2 and decrease k by 1;
        ## and [W3] cut all edges in E_2 in F_2 and decrease k by |E_2|
        pass
    def branching_rule_4(self, X):
        """Case 4. All leaves in X_2 are in the same tree in F_2 but no two are siblings"""
        ## Assumption: there are no sibling sets.
        ## Let u_2, v_2 be a elements of X_2 such that the parent p_2 of u_2 has degree 2 in F_2[l(X_2)].
        ## (l(X_2) is the set of labels for leaves X_2)
        ## (F_2[l(X_2)] is presumably the tree prune of all leaves except X_2; subforest of F_2 induced by l(X_2))
        for i in range(len(X)-1):
            u = X[i]
            for j in range(i+1,len(X)):
                v = X[j]
                if self.is_sibling_set([u, v]):
                    return 0 ## sibling sets violate the premise
        #### TODO: validate that the bit of code above is in line with the algorithm
        ## We split this case into three subcases depending on the size |X_1|
        ## and the number of unlabelled vertices on a path connecting two leaves in X_2.
        ## Let u_2 and v_2 be two arbitrary leaves in X_2, and let P be the path in F_2 that connects u_2 and v_2.
        ## Denote by int(P) = {w_1,...,w_h}, h>=2, the set of vertices on P that are unlabelled.
        ## Subcase 4.1: the path P consists of at least five vertices, i.e. h>=3. (Fig. 2(B))
        ## Branching rule 4.1: letu_2 and v_2 be leaves in X_2 such that the path P = {u_2,w_1,...,w_h,v_2} in F_2
        ## satisfies h>=3. Then branch into (h+2) eays:
        ## [W1] cut u_2 in F_2 and decrease k by 1;
        ## [W2] but v_2 in F_2 and decrease k by 1;
        ## and, for each 1<=i<=h, [W(2+1)] cut the edges incident to int(P)\{w_i} but not on the path P, and decrease k
        ## by the number of edges cut.
        ## Subcase 4.2: E_2 = {e_1,...,e_h} with e>=3. (Fig. 2(C))
        ## Brancing rule 4.2: branch into (2+h) ways:
        ## [W1] cut u_2 in F_2 and decrease k by 1;
        ## [W2] cut v_2 in F_2 and decrease k by 1;
        ## and, for each 1<=i<=h, [W(2+i)] for the edge e_i in E_2, cut all edges in E_2 except e_i and decrease k by h-1
        ## Subcase 4.3: E_2={e_1,e_2}. (Fig. 2(D))
        ## Branching rule 4.3: decrease k by 1, and branch into three ways:
        ## [W0] cut u_2 in F_2;
        ## [W1] cut e_1 in E_2;
        ## and [W2] cut e_2 in E_2.
        pass
    # def E_double_prime(self, X):
    #     """
    #     Sibling set X has parent p and the degree of p is at least |X|+2.
    #     Let E" be the set of edges that are incident to the parent p of X but not incident to the leaves in X.
    #     """
    #     if not self.is_sibling_set(X): return []
    #     parent = self.get_parent(X) ## 1 parent
    #     X_indices = set(int(v) for v in X)
    #     output = []
    #     for e in parent.all_edges():
    #         ## skip if edge is connected to any vertex in X
    #         if int(e.source()) in X_indices or int(e.target()) in X_indices: continue
    #         output.append(e)
    #     return output
    # def meta_step_2(self, X):
    #     """Case 2. X_2 is a sibling set in F_2."""
    #     if (not self.is_bss(X)
    #         and self.is_sibling_set(X)):
    #         E_double_prime = self.E_double_prime(X)
    #     else:
    #         return 0 ## operation not executed
    #     ## get edges between X and parent
    #     X_indices = set(int(v) for v in X)
    #     parent = self.get_parent(X)
    #     E_X_parent = [e for e in parent.all_edges() if
    #                   (int(e.source()) in X_indices
    #                    or int(e.target()) in X_indices)]
    #     ## meta-step 2.1
    #     ## if |E"| > |X_2|, then pick any set E_1 of |X_2|-1 edges incident to the leaves in X_2,
    #     ## and pick any set E_2 of |X_2| edges in E", remove all edges in E_1 union E_2
    #     if len(E_double_prime) > len(X):
    #         E1 = E_X_parent[:-1] ## pick any set E_1 of |X_2|-1 edges incident to leaves in X
    #         E2 = E_double_prime[:len(X)] ## pick any set E_2 of |X_2| edges in E"
    #         for e in E1 + E2:
    #             self.remove_edge(e)
    #     ## meta-step 2.2
    #     ## if |E"| <= |X_2|, then pick any set E'_1 of |E"|-1 edges incident to the leaves in X_2,
    #     ## and remove all edges in E'_1 union E".
    #     else:
    #         E_prime_1 = E_X_parent[:(len(E_double_prime)-1)] ## pick any set E'_1 of |E"|-1 edges incident to the leaves in X_2
    #         for e in E_prime_1 + E_double_prime:
    #             self.remove_edge(e)
    #     return 1 ## operation executed successfully
    # @temporarily_undirected
    # def meta_step_3(self, X, preserve_direction = False):
    #     """Case 3. X_2 contains a sibling set Y_2 with |Y_2|>=2. Assumption: X_2 is not a sibling set."""
    #     ## Assumption: X is not a sibling set. assume also: graph is acyclic
    #     ## Find a sibling set
    #     ## 'Let u_2, v_2 be elements of Y_2'
    #     if len(X) < 3: return 0 ## operation not executed. Too few vertices
    #     Y = None
    #     break_loop = False
    #     for i in range(len(X)-1):
    #         u = X[i]
    #         for j in range(i+1,len(X)):
    #             v = X[j]
    #             if self.is_sibling_set([u, v]):
    #                 Y = [u, v]
    #                 break_loop = True
    #                 break
    #         if break_loop: break
    #     if Y is None: return 0 ## operation not executed. Premise that X_2 contains a sibling set Y_2 is false
    #     ## 'Let w_2 be a leaf in X_2 that is not a sibling of u_2 and v_2.'
    #     Y_parent = self.get_parent(Y) ## 1 parent
    #     for w in X:
    #         w_neighbour_indices = [int(n) for n in w.all_neighbours()]
    #         if Y_parent not in w_neighbour_indices: break ## get first vertex without u+v's parent as a neighbour
    #     ## 'Let e be an edge incident to the parent of w_2 but not on the path between u_w and w_2.'
    #     # g_is_directed = self.is_directed()
    #     # self.set_directed(False) ## set undirected to simplify steps
    #     path_vertices = gt.shortest_path(self, u, w)[0]
    #     path_vertices_index = [int(n) for n in path_vertices]
    #     w_parent = path_vertices[-2] ## -1 is w
    #     for n in w_parent.all_neighbours():
    #         ## check if n is along path between u & w
    #         if int(n) not in path_vertices_index:
    #             ## if not along path, get edge between n & w_parent
    #             e = self.edge(n, w_parent)
    #             break
    #     ## 'Then remove the edge e, the edge e_u incident to u, and the edge e_w incident to w'
    #     self.remove_edge(e)
    #     self.remove_edge(self.edge(u, Y_parent))
    #     self.remove_edge(self.edge(w, w_parent))
    #     # ## set direction back to original (even though this operation only really works on undirected graphs anyway...)
    #     # self.set_directed(g_is_directed)
    #     return 1 ## operation executed successfully
    # @temporarily_undirected
    # def meta_step_4(self, X, preserve_direction = False):
    #     """Case 4. All leaves in X_2 are in the same tree in F_2 but no two are siblings"""
    #     ## Assumption: there are no sibling sets.
    #     ## Let u_2, v_2 be a elements of X_2 such that the parent p_2 of u_2 has degree 2 in F_2[l(X_2)].
    #     ## (l(X_2) is the set of labels for leaves X_2)
    #     ## (F_2[l(X_2)] is presumably the tree prune of all leaves except X_2; subforest of F_2 induced by l(X_2))
    #     for i in range(len(X)-1):
    #         u = X[i]
    #         for j in range(i+1,len(X)):
    #             v = X[j]
    #             if self.is_sibling_set([u, v]):
    #                 return 0 ## sibling sets violate the premise
    #     ## Let e be an edge in F_2 that is incident to p_2 but not on the path P_{uv} between u_2 and v_2.
    #     # g_is_directed = self.is_directed()
    #     # self.set_directed(False) ## set undirected to simplify steps
    #     path_vertices = gt.shortest_path(self, u, v)[0]
    #     path_vertices_index = [int(n) for n in path_vertices]
    #     p = path_vertices[1] ## 0 is u, -1 is v
    #     e = None
    #     for n in p.all_neighbours():
    #         ## check if n is along path between u & w
    #         if int(n) not in path_vertices_index:
    #             ## if not along path, get edge between n & p
    #             e = self.edge(n, p)
    #             break
    #     ## Then cut the edge e, the edge e_u incident to u_2, and the edge e_v incident to v_2.
    #     if e is None: return 0
    #     # try:
    #     #     self.remove_edge(e)
    #     # except Exception as er:
    #     #     print(path_vertices)
    #     #     print(p.all_neighbours())
    #     #     print([int(n) for n in path_vertices])
    #     #     print(p)
    #     #     print(n)
    #     #     raise er
    #     self.remove_edge(self.edge(u, p))
    #     self.remove_edge(self.edge(v, path_vertices[-2]))
    #     # ## set direction back to original (even though this operation only really works on undirected graphs anyway...)
    #     # self.set_directed(g_is_directed)
    #     return 1 ## operation executed successfully
    # def meta_step_switch(self, X, preserve_direction = False):
    #     """Switch for meta-steps"""
    #     executed = self.meta_step_1(X) ## if leaves are not in same tree
    #     if not executed: executed = self.meta_step_2(X) ## if X_2 is a sibling set
    #     if not executed: executed = self.meta_step_3(X) ## if X_2 contains a sibling set
    #     if not executed: executed = self.meta_step_4(X) ## if all leaves are in same tree but no 2 are siblings
    #     ## return whether any meta-steps were executed (0/1)
    #     return executed
