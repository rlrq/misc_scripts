#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines
from gff_manip import GFF

def extract_features_and_subfeatures(gff_fname, feature_id_fname, fout, fout_fmt = "GFF",
                                     fmt = None, attr_mod = '', memsave = False):
    print("Reading files")
    if fmt is None:
        fmt = "BED" if gff_fname[-3:].upper() == "BED" else "GFF"
    gff = GFF(fname = gff_fname, fmt = fmt, attr_mod = attr_mod, memsave = memsave)
    feature_ids = splitlines(feature_id_fname)
    print("Retrieving features")
    output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
    print("Writing features")
    gff.write_i(fout, output_features, fmt = fout_fmt)
    return

parser = argparse.ArgumentParser(description="subset GFF file")
parser.add_argument("gff", type=str, help="path to GFF3/BED file")
parser.add_argument("feature_id_fname", type=str, help="path to file containing feature IDs")
parser.add_argument("fout", type=str, help="output file")
parser.add_argument("fout_fmt", type=str, help="output format [GFF/BED]", nargs='?', default="GFF")
parser.add_argument("--memsave", help="memory-saving mode", action="store_true", default=False)
parser.add_argument("--fmt", type=str, help="input format [GFF/BED]")
parser.add_argument("--attr-mod", type=str, help="attribute modification")
args = parser.parse_args()

extract_features_and_subfeatures(args.gff, args.feature_id_fname, args.fout, fout_fmt = args.fout_fmt,
                                 memsave = args.memsave, fmt = args.fmt, attr_mod = args.attr_mod)

# gff_fname, feature_id_fname, fout = sys.argv[1:4]
# fout_fmt = "GFF" if len(sys.argv) < 5 else sys.argv[4]

# extract_features_and_subfeatures(gff_fname, feature_id_fname, fout, fout_fmt = fout_fmt)
