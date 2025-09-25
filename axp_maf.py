import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import newick

## based on https://www.sciencedirect.com/science/article/pii/S0304397514008238#br0060
## Parameterized and approximation algorithms for maximum agreement forest in multifurcating trees
## 2015; Jianer Chen, Jia-Hao Fan, & Sing-Hoi Sze; Theoretical Computer Science. 562:496-512

#### for development, use these trees, which have some multifurcating subtrees
from Bio import SeqIO
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
    '''
    t is a tree object read by package newick.
    a is the attribute name in which to store the output of f.
    f is a function that takes a node and outputs a value.
    '''
    nodes = t.get_leaves()
    for n in nodes:
        setattr(n, a, f(n))
    return

def set_leaf_attributes(t, a_f_dict):
    '''
    t is a tree object read by package newick.
    a_f_dict is a dictionary of {"<attribute name>": f}, where
      f is a function that takes a node and outputs a value to
      be stored in attribute with name <attribute name>.
    '''
    for a, f in a_f_dict.items():
        set_leaf_attribute(t, a, f)
    return

## get seqid of longest transcript for each gene
seqid_gid = {seqid: re.search("(?<=\\|).+(?=\\|CDS\\|)", seqid).group(0) for seqid in seqlens}
gid_seqid = basics.invert_dict(seqid_gid)
seqid_longest = {gid: max(seqids, key = lambda seqid:seqlens[seqid]) for gid, seqids in gid_seqid.items()}
## get longest gene for each accession (per nwk)
seqid_acc = {seqid: ("TAIR10" if "Col-0_ref" in seqid
                     else re.search("(?<=Reference\\|)[^|]+", seqid).group(0) if re.search("\\d\\|G\\d", seqid)
                     else re.search("(?<=Reference\\|).+(?=\\.g[-\\d]+?\\|CDS)", seqid).group(0))
             for seqid in seqlens}
longest_per_acc_in_tree = lambda t:list({acc: max(seqids, key = lambda seqid:seqlens[seqid]) for acc, seqids in basics.invert_dict({leaf_name: seqid_acc[leaf_name] for leaf_name in t.get_leaf_names()}).items()}.values())

## read trees
t1 = newick.read(nwk1)[0]
t1.prune_by_names(longest_per_acc_in_tree(t1), inverse = True)
t1.remove_redundant_nodes(keep_leaf_name = True)
t2 = newick.read(nwk2)[0]
t2.prune_by_names(longest_per_acc_in_tree(t2), inverse = True)
t2.remove_redundant_nodes(keep_leaf_name = True)
## remove nodes from accessions that aren't shared between the two trees

## new functions, replace hardcoded dicts above
get_gid = lambda seqid: re.search("(?<=\\|).+(?=\\|CDS\\|)", seqid).group(0)
get_accid = lambda seqid: ("TAIR10" if "Col-0_ref" in seqid
                           else re.search("(?<=Reference\\|)[^|]+", seqid).group(0) if re.search("\\d\\|G\\d", seqid)
                           else re.search("(?<=Reference\\|).+(?=\\.g[-\\d]+?\\|CDS)", seqid).group(0))
attr_map = {"seqlen": lambda n:seqlens[n.name],
            "accid": lambda n:seqid_acc[n.name],
            "gid": lambda n:seqid_gid[n.name]}

## read trees v2
t1 = newick.read(nwk1)[0]
set_leaf_attributes(t1, attr_map)
t2 = newick.read(nwk2)[0]
set_leaf_attributes(t2, attr_map)



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

def reduction_rule_1(t):
    '''
      If a label l is in a single-vertex tree in one of the forests F_1 and F_2,
    then remove the edge (if any) incident to the label l in the other forest.
    '''
    pass

def reduction_rule_2(t):
    '''
      Let X_1 be a BSS in F_1, and let X_2 be a set of leaves in F_2 such that
    l(X_1) = l(X_2). If X_2 is the leaf set of a single-edge tree in F_2, or if
    X_2 is a sibling set whose parent has degree at most |X_2|+1, then shrink
    X_1 in F_1, and shrink X_2 in F_2.
      Note that after applying Reduction Rule 2, if a vertex v in either F_1 or
    F_2 is adjacent to the new leaf (there is at most one such vertex), then the
    degree of v is at most 2. In particular, if v is unlabelled, then v will be
    contracted.
    '''
    pass

