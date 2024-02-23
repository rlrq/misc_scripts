#!/usr/bin/python3

import sys
import copy
import numpy
import logging
import pathlib
import argparse
import dendropy

# sys.path.append("/mnt/chaelab/rachelle/src")
# from basics import head, tail

logging.basicConfig(level=logging.INFO)


###########################
##    FUNCTIONS COPIED   ##
##  FROM split_by_decay  ##
###########################

def edge_len(e):
    return (e.length if e.length is not None else 0)

def node_label(n):
    return n.taxon.label

def branch_lengths(t):
    ## parse root edges (since we're assuming unrooted tree)
    root_edges = tuple(t.seed_node.child_edges())
    num_root_edges = len(root_edges)
    root_edges_lengths = ([0] if num_root_edges < 2 else
                          [sum([edge_len(e) for e in root_edges])] if num_root_edges == 2 else
                          [edge_len(e) for e in root_edges])
    ## parse other edges
    branch_lengths = [edge_len(e) for e in t.preorder_edge_iter() if
                      (e is not t.seed_node.edge and not e in root_edges)] + root_edges_lengths
    return branch_lengths

def max_branch_length(t):
    return max(branch_lengths(t))

def sd_mean(t):
    ## parse root edges (since we're assuming unrooted tree)
    root_edges = tuple(t.seed_node.child_edges())
    if len(root_edges) == 2:
        root_edges_lengths = [sum(edge_len(e) for e in root_edges)]
    else:
        root_edges_lengths = [edge_len(e) for e in root_edges]
    # root_edges_length = (0 if len(root_edges) != 2 else sum([edge_len(e) for e in root_edges]))
    ## parse other edges
    branch_lengths = [edge_len(e) for e in t.preorder_edge_iter() if
                      (e is not t.seed_node.edge and not e in root_edges)] + root_edges_lengths
    ## return NaN if only 1 leaf or all branch lengths == 0
    if (len(branch_lengths) <= 1 or all(map(lambda l: l==0, branch_lengths))):
        return float("NaN")
    sdmean = numpy.std(branch_lengths)/numpy.mean(branch_lengths)
    logging.debug(f"{sdmean} {branch_lengths}")
    return sdmean

def parse_tree(tree, rooted = False, schema = None):
    ## figure out schema
    if schema is None:
        ext = tree.split('.')[-1].lower()
        schema = ("newick" if ext in {"newick", "nwk"} \
                  else "nexus" if ext in {"nexus", "nex", "nxs"} \
                  else "nexml" if ext in {"nexml"} \
                  else "phylip" if ext in {"phylip", "phy"} \
                  else "fasta" if ext[0] == 'f' \
                  else "newick") ## apprently newick also uses '.tree' extension?
    ## read tree
    t = (dendropy.Tree.get(path = tree, schema = schema, rooting = "force-rooted") if rooted
         else dendropy.Tree.get(path = tree, schema = schema, rooting = "default-unrooted"))
    return t

def split_at_longest_edge(t, min_size: int = 1):
    if min_size < 1: min_size = 1
    ## some functions
    edge_child = lambda e: e.head_node
    ## get edges with length (exclude root)
    edges = [e for e in t.preorder_edge_iter() if e is not t.seed_node.edge]
    ## if edges have no length or is terminal leaf, return tree
    if (sum(branch_lengths(t)) == 0 or \
        len(t.leaf_nodes()) <= min_size):
        print(branch_lengths(t))
        return [t, None]
    ## parse root edges
    root_edges = tuple(t.seed_node.child_edges())
    root_edges_length = (0 if len(root_edges) != 2 else sum([edge_len(e) for e in root_edges]))
    ## get longest edge
    longest_edge = max(edges, key = lambda e: edge_len(e))
    ## instantiate trees to be returned
    t1 = copy.deepcopy(t)
    t2 = copy.deepcopy(t)
    ## if sum of root edges are longer than longest_edge, split at root
    if root_edges_length >= edge_len(longest_edge):
        t1.retain_taxa_with_labels([node_label(n) for n in edge_child(root_edges[0]).leaf_nodes()])
        t2.retain_taxa_with_labels([node_label(n) for n in edge_child(root_edges[1]).leaf_nodes()])
    else:
        t1_taxa = [node_label(n) for n in edge_child(longest_edge).leaf_nodes()]
        t1.retain_taxa_with_labels(t1_taxa)
        t2.retain_taxa_with_labels([node_label(n) for n in t.leaf_nodes() if node_label(n) not in t1_taxa])
    return [t1, t2]


######################
##  MAIN FUNCTIONS  ##
######################

## members should be comma-separated leaf names
def report_separation(fout, tree, og, grp1, grp2, members = '', schema = None):
    t = parse_tree(tree, schema = schema, rooted = False)
    longest_edge = max_branch_length(t)
    sdmean = sd_mean(t)
    t1, t2 = split_at_longest_edge(t, min_size = 1)
    members = set(members.split(','))
    t1_leaves = set(node_label(n) for n in t1.leaf_nodes())
    t2_leaves = set(node_label(n) for n in t2.leaf_nodes())
    status = 'T' if (members == t1_leaves or members == t2_leaves) else 'F'
    sdmean_t1 = sd_mean(t1)
    sdmean_t2 = sd_mean(t2)
    with open(fout, "a+") as f:
        logging.debug(fout)
        ## columns: og, grp1, grp2, whether longest edge separates groups, length of longest edge, sd/mean
        f.write(f"{og}\t{grp1}\t{grp2}\t{status}\t{longest_edge}\t{sdmean}\n")
        # ## columns: og, grp1, grp2, whether longest edge separates groups, length of longest edge, sd/mean, sd/mean of tree1, sd/mean of tree2
        # f.write(f"{og}\t{grp1}\t{grp2}\t{status}\t{longest_edge}\t{sdmean}\t{sdmean_t1}\t{sdmean_t2}\n")
    return

## create parser
parser = argparse.ArgumentParser(description="split leaves of tree into subgroups based on decay")

## add arguments
parser.add_argument("-o", "--out", type = pathlib.Path, help = "path to output file")
parser.add_argument("-t", "--tree", type = pathlib.Path, help = "path to phylogenetic tree")
parser.add_argument("--schema", type = str, default = None, help = "tree schema (e.g. nexus, newick)")
parser.add_argument("--og", type = str, help = "og name")
parser.add_argument("--g1", "--grp1", type = str, dest = "grp1", help = "group 1 id (integer value)")
parser.add_argument("--g2", "--grp2", type = str, dest = "grp2", help = "group 2 id (integer value)")
parser.add_argument("-m", "--members", "--grp1members", "--grp2members",
                    type = str, dest = "members", help = "comma-separated leaf names for either group")

args = parser.parse_args()
if len(sys.argv) > 1:
    report_separation(args.out, tree = args.tree, schema = args.schema,
                      og = args.og, grp1 = args.grp1, grp2 = args.grp2,
                      members = args.members)
