#!/usr/bin/python3
import argparse

### TODO

# import sys
# sys.path.append("/mnt/chaelab/rchelle/src")

# import annotation as ann

# def read_annotation(fname, fmt = None, attr_mod = '', memsave = False):
#     print("Reading annotations")
#     if fmt is None:
#         fmt = "BED" if fname[-3:].upper() == "BED" else "GFF"
#     gff = ann.GFF(fname = fname, fmt = fmt, attr_mod = attr_mod, memsave = memsave)
#     return gff



# parser = argparse.ArgumentParser(description="identify head-to-head genes within a given gene set")
# parser.add_argument("gff", type=str, help="path to GFF3/BED file")
# parser.add_argument("feature_id_fname", type=str, help="path to file containing feature IDs")
# parser.add_argument("fout", type=str, help="output file")
# parser.add_argument("--memsave", help="memory-saving mode", action="store_true", default=False)
# parser.add_argument("--fmt", type=str, help="input format [GFF/BED]")
# parser.add_argument("--attr-mod", type=str, help="attribute modification")
# # ## if self is true, then only the features themselves (no parents, no children) will be returned
# # ## if both parent and child are false, then full features (feature + parents + children) will be returned
# # ## if both parent and child are false and include_cousins is true, then subfeatures related through parent features will be included
# # group = parser.add_mutually_exclusive_group()
# # group.add_argument("--parent", help="get feature and parent features only", action="store_true", default=False)
# # group.add_argument("--child", help="get feature and subfeatures only", action="store_true", default=False)
# # group.add_argument("--self", help="get feature only", action="store_true", default=False)
# # group.add_argument("--include-cousins", help="get feature and parent features and all subfeatures of parent features",
# #                    action="store_true", default=False)
# args = parser.parse_args()
