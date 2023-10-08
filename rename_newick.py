#!/usr/bin/python3
import argparse
from Bio import Phylo

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines

def rename_nodes(fin, fout, names):
    ## parse name mapping
    names = [x.split('\t') for x in splitlines(names)]
    names = {x[0]: x[1] for x in names}
    ## parse tree
    t = Phylo.read(fin, "newick")
    ## get all nodes and leaves
    nodes = t.get_nonterminals() + t.get_terminals()
    for node in nodes:
        if node.name in names:
            node.name = names[node.name]
    ## write
    Phylo.write(t, fout, "newick")

parser = argparse.ArgumentParser(description="rename leaves of newick tree")
parser.add_argument("nwk", type=str, help="path to input newick file")
parser.add_argument("names", type=str, help="path to tab separated file of [old name] [new name]")
parser.add_argument("fout", type=str, help="path to output newick file", nargs='?', default=None)
parser.add_argument("-i", help="in-place", action="store_true", default=False)
args = parser.parse_args()

fout = (args.fout if not args.i else args.nwk)
rename_nodes(args.nwk, fout, args.names)
