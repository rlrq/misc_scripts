#!/usr/bin/python3

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines
from gff_manip import GFF

def extract_features_and_subfeatures(gff_fname, feature_id_fname, fout, fout_fmt = "GFF", fmt = "GFF"):
    print("Reading files")
    if gff_fname[-3:].upper() == "BED": fmt = "BED"
    gff = GFF(fname = gff_fname, fmt = fmt)
    feature_ids = splitlines(feature_id_fname)
    print("Retrieving features")
    output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
    print("Writing features")
    gff.write_i(fout, output_features, fmt = fout_fmt)
    return

gff_fname, feature_id_fname, fout = sys.argv[1:4]
fout_fmt = "GFF" if len(sys.argv) < 5 else sys.argv[4]

extract_features_and_subfeatures(gff_fname, feature_id_fname, fout, fout_fmt = fout_fmt)
