#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines
from gff_manip import GFF

def extract_features_and_subfeatures(gff_fname, feature_id_fname, fout, fout_fmt = "GFF",
                                     parent = False, child = False, self_only = False,
                                     include_cousins = False,
                                     fmt = None, attr_mod = '', memsave = False):
    print("Reading files")
    if fmt is None:
        fmt = "BED" if gff_fname[-3:].upper() == "BED" else "GFF"
    gff = GFF(fname = gff_fname, fmt = fmt, attr_mod = attr_mod, memsave = memsave)
    feature_ids = splitlines(feature_id_fname)
    print(feature_ids)
    print("Retrieving features")
    if self_only:
        output_features = gff.get_id(feature_ids, index = True, output_list = True)
    elif parent:
        output_features = gff.get_features_and_parent_features(feature_ids, index = True, full = True)
    elif child:
        output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
    else:
        output_features = gff.get_related_features(feature_ids, index = True, full = True, cousins = include_cousins)
    output_features = sorted(set(output_features)) ## sort and remove repeat indices
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
## if self is true, then only the features themselves (no parents, no children) will be returned
## if both parent and child are false, then full features (feature + parents + children) will be returned
## if both parent and child are false and include_cousins is true, then subfeatures related through parent features will be included
group = parser.add_mutually_exclusive_group()
group.add_argument("--parent", help="get feature and parent features only", action="store_true", default=False)
group.add_argument("--child", help="get feature and subfeatures only", action="store_true", default=False)
group.add_argument("--self", help="get feature only", action="store_true", default=False)
group.add_argument("--include-cousins", help="get feature and parent features and all subfeatures of parent features",
                   action="store_true", default=False)
args = parser.parse_args()

extract_features_and_subfeatures(args.gff, args.feature_id_fname, args.fout, fout_fmt = args.fout_fmt,
                                 parent = args.parent, child = args.child, self_only = args.self,
                                 include_cousins = args.include_cousins,
                                 memsave = args.memsave, fmt = args.fmt, attr_mod = args.attr_mod)

# gff_fname, feature_id_fname, fout = sys.argv[1:4]
# fout_fmt = "GFF" if len(sys.argv) < 5 else sys.argv[4]

# extract_features_and_subfeatures(gff_fname, feature_id_fname, fout, fout_fmt = fout_fmt)
