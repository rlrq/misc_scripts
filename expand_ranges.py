#!/usr/bin/python3

"""
Expands ranges in GFF3 or BED file.
Bound by [0, inf) if BED fmt or [1, inf) if GFF fmt.
"""

import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines

def expand_range(fname, fout, distance = 0, fmt = None):
    print("Reading files")
    if fmt is None:
        fmt = "BED" if fname[-3:].upper() == "BED" else "GFF"
    else:
        fmt = fmt.upper()
    lower_limit = 0 if fmt == "BED" else 1
    start_i = 1 if fmt == "BED" else 3
    end_i = 2 if fmt == "BED" else 4
    dat = [x.split('\t') for x in splitlines(fname)]
    print("Expanding ranges")
    for entry in dat:
        if entry[0][0] == '#':
            if entry == "##FASTA": break
            else: continue
        entry[start_i] = max(lower_limit, int(entry[start_i]) - distance)
        entry[end_i] = int(entry[end_i]) + distance
    print("Writing features")
    with open(fout, "w+") as f:
        for entry in dat:
            silence = f.write('\t'.join(map(str, entry)) + '\n')
    return

parser = argparse.ArgumentParser(description="expand ranges in BED/GFF file")
parser.add_argument("fname", type=str, help="path to BED/GFF3 file")
parser.add_argument("fout", type=str, help="output file")
parser.add_argument("distance", type=int, help="distance to expand (bp)")
# parser.add_argument("fout_fmt", type=str, help="output format [GFF/BED]", nargs='?', default="GFF")
# parser.add_argument("--memsave", help="memory-saving mode", action="store_true", default=False)
parser.add_argument("--fmt", type=str, help="input format [GFF/BED]")
# parser.add_argument("--attr-mod", type=str, help="attribute modification")
args = parser.parse_args()

expand_range(args.fname, args.fout, distance = args.distance, fmt = args.fmt)
