#!/usr/bin/python3
import argparse
from Bio import Phylo

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines

def collapse_groups(fin, fout, fgroups):
    ## parse file mapping leaves to groups
    mapping = [x.split('\t') for x in splitlines(fgroups)]
    print(mapping)
    mapping = {x[0]: x[1] for x in mapping}
    def map_label(s):
        return mapping.get(s, s)
    ## parse tree
    t = Phylo.read(fin, "newick")
    ## iterate through all sub-trees
    def traverse(t):
        if t.is_terminal():
            ## rename
            t.name = map_labels(t.name)
            return {t.name}
        else:
            ## sort out group of sub trees
            groups = set()
            for child in t:
                groups.union(traverse(child))
            ## if all descendant leaves are in same group,
            ## prune leaves and make this node a leaf instead
            if len(groups) == 1:
                leaves = t.get_terminals()
                t.name = leaves[0].name
                for leaf in leaves:
                    t.prune(leaf)
            return
    ## write
    Phylo.write(t, fout, "newick")
    return

# def collapse_clade(fin, fout, features):
#     ## parse file mapping leaves to groups
#     mapping = [x.split('\t') for x in splitlines(features)]
#     mapping = {x[0]: x[1] for x in mapping}
#     def map_label(s):
#         return mapping.get(s, s)
#     ## get index of end of node data
#     def end_of_node(s):
#         ## get index of end of tip data
#         i_end = -1
#         try:
#             i_end = s.index(',')
#         except ValueError:
#             i_end = s.index(')')
#         return i_end
#     ## get tip name
#     def tip_name(s):
#         ## get tip data
#         s_tip = s[:end_of_node(s)]
#         ## get end of tip name
#         i = len(s_tip)-1
#         while s_tip[i] != ':' and i >= 0:
#             i -= 1
#         ## get tip name
#         if i < 0: return s_tip
#         else: return s_tip[:i]
#     ## process a tree
#     def condense(s):
#         leaves = s.split(',')
#         labels = [tip_name(leaf) for leaf in leaves]
#         mapped = [map_label(label) for label in labels]
#         groups = set(mapped)
#         if len(groups) == 1:
#             return mapped[0], groups
#         else:
#             return ','.join([leaf.replace(labels[i], mapped[i]) for i, leaf in enumerate(leaves)]), groups
#     ## get subtree (without parentheses) (given index of closing parenthesis)
#     def get_subtree(s, i):
#         start = i-1
#         while start > 0 and s[start] != '(':
#             start -= 1
#         return s[start+1:i]
#     ## recursively map tip names
#     def traverse(s):
#         i = 0
#         while i < len(s):
#             ## if end of subtree, extract subtree
#             if s[i] == ')':
#                 subtree = get_subtree(s, i)
#                 nwk, groups = condense(subtree)
#                 s = s[:i-len(subtree)] + nwk + s[i:]
#         ## if the first node points to a subtree, recur
#         elif s[0] == '(':
#             node_children, node_nwk = helper(s[1:])
#             return node_children, '(' + node_nwk
#         ## if the first node points to a leaf
#         else:
#             ## tip1
#             eon = end_of_node(s) ## end of node
#             label_1 = tip_name(s)
#             node_meta = s[len(label_1):]
#             new_label = map_label(label_1)
#             curr_node_children = {new_label}
#             curr_node_nwk = new_label + node_meta + s[eon]
#             ## subsequent tips
#             next_node_children, next_node_nwk = helper(s[end_of_node(s)+1:])
#             if len(curr_node_childre.union(next_node_children)) == 1:
#                 return curr_node_children, new_label + ":0.0"

parser = argparse.ArgumentParser(description="collapse clades of a newick tree where all descendants belong to the same group/phylum and rename all leaves by group/phylum name")
parser.add_argument("nwk", type=str, help="path to input newick file")
parser.add_argument("groups", type=str, help="path to tab separated file of [leaf name] [group name]")
parser.add_argument("fout", type=str, help="path to output newick file", nargs='?', default=None)
parser.add_argument("-i", help="in-place", action="store_true", default=False)
args = parser.parse_args()

fout = (args.fout if not args.i else args.nwk)
collapse_groups(args.nwk, fout, args.groups)
