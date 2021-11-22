#!/usr/bin/python3

import sys
import copy
import numpy
import pathlib
import argparse
import dendropy
import logging

logging.basicConfig(level=logging.WARNING)

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
    root_edges_length = (0 if len(root_edges) != 2 else sum([edge_len(e) for e in root_edges]))
    ## parse other edges
    branch_lengths = [edge_len(e) for e in t.preorder_edge_iter() if
                      (e is not t.seed_node.edge and not e in root_edges)] + [root_edges_length]
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
    if (sum([edge_len(e) for e in edges]) == 0 or \
        len(t.leaf_nodes()) <= min_size):
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

def split_by_decay(tree, ratio = 2, rooted = False, min_size: int = 1, return_tree = False,
                   return_taxa = False, min_split_len = 0.02, log_ratio = 2, max_keep_len = 0.1,
                   # last_ratio = None, gradient = -0.1, last_gradient = None,
                   **kwargs):
    """
    ratio: max sd/mean threshold allowed for clade
    """
    def sd_mean_terminal(t):
        num_branches = len(tuple(t.preorder_edge_iter()))
        threshold = max(ratio, log_ratio/numpy.log10(num_branches)) ## this is arbitrary
        return sd_mean(t) <= threshold
    def split_recursively(t):
        logging.debug("iter")
        def fmt_return(t):
            if return_tree: return [t] ## return tree obj
            elif return_taxa: return [[node_label(n) for n in t.leaf_nodes()]] ## return list of taxa in tree
            else: return [[sd_mean(t), [node_label(n) for n in t.leaf_nodes()]]] ## sd/mean & taxa
        if t is None: return []
        else:
            max_branch_len = max_branch_length(t)
            if (max_branch_len <= max_keep_len and
                (len(t.leaf_nodes()) <= min_size or
                 max_branch_len < min_split_len or
                 sd_mean_terminal(t))):
                return fmt_return(t)
            else:
                output = []
                ## start splitting
                t1, t2 = split_at_longest_edge(t, min_size = min_size)
                if t2 is None:
                    output.extend(fmt_return(t1))
                else:
                    output.extend(split_recursively(t1))
                    output.extend(split_recursively(t2))
                return output
    t = parse_tree(tree, rooted = rooted, **kwargs)
    ## start splitting
    ## ...actually I don't think rooted vs unrooted trees make a difference lol
    ## I don't think there's a need to re-root at all actually.
    return split_recursively(t)

def split_by_decay_and_write(fout, **kwargs):
    groups = split_by_decay(return_tree = False, return_taxa = False, **kwargs)
    with open(fout, "w+") as f:
        logging.info(fout)
        for i, dat in enumerate(groups):
            ## write <grp number>\t{sd/mean}\t{grp members}
            f.write(f"{i+1}\t{dat[0]}\t{','.join(dat[1])}\n")
    return


## create parser
parser = argparse.ArgumentParser(description="split leaves of tree into subgroups based on decay")

## add arguments
parser.add_argument("-t", "--tree", type = pathlib.Path, help = "path to phylogenetic tree")
parser.add_argument("--schema", type = str, default = None, help = "tree schema (e.g. nexus, newick)")
parser.add_argument("-r", "--ratio", type = float, default = 2,
                    help = "sd/mean threshold for clade definition")
parser.add_argument("-l", "--logratio", type = float, default = 3,
                    help = "1/log10 threshold for clade definition")
parser.add_argument("-o", "--out", type = pathlib.Path, help = "path to output file")
parser.add_argument("--minsize", type = int, help = "minimum clade size", default = 1)
parser.add_argument("--minsplitlen", type = float, help = "minimum eligible length of branch for splitting", default = 0.02) ## if a tree does not contain branches of this length or greater, the tree will not be split further
parser.add_argument("--maxkeeplen", type = float, help = "maximum eligible length of branch to be conserved; branches longer than this will be split regardless of sd/mean", default = 0.1) ## if a tree does not contain branches of this length or greater, the tree will not be split further
# parser.add_argument("-g", "--gradient", type = float,
#                     help = "minimum non-positive decay gradient for clade definition")
# parser.add_argument("--rooted", action = "store_true", help = "input is rooted tree", default = False)

args = parser.parse_args()
# print(args)
split_by_decay_and_write(args.out, tree = args.tree, schema = args.schema,
                         ratio = args.ratio, log_ratio = args.logratio,
                         min_size = args.minsize, min_split_len = args.minsplitlen,
                         max_keep_len = args.maxkeeplen)
            
# tree_p = "/mnt/chaelab/rachelle/pav_haplotype/results/og_decay/og153_1/og153_1.CDS.mafft.nwk"
# t = dendropy.Tree.get(path = "/mnt/chaelab/rachelle/pav_haplotype/results/og_decay/og153_1/og153_1.CDS.mafft.nwk", schema = "newick")
# x = split_by_decay(tree_p)

# tree_str = "(dog:20, (elephant:30, horse:60):20);"
# t_d = dendropy.Tree.get_from_string(tree_str, schema = "newick", rooting = "default-unrooted")
# [x.length for x in t_d.preorder_edge_iter()]
# print(t_d.as_ascii_plot())
# edge_d = min(t_d.preorder_edge_iter().key = lambda x: (x.length if x.length else float("Inf")))
# t_d.reroot_at_edge(edge_d, update_bipartitions = True)
# print(t_d.as_ascii_plot())

# import io
# from Bio import Phylo
# t_b = Phylo.read(io.StringIO(tree_str), "newick")
# Phylo.draw_ascii(t_b)
