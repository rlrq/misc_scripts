#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from gff_manip import GFF

## parses gff files and writes tsv file of <child>\t<parent>
## output will be sorted
def get_and_write_gff_relationships(gff_fname, fout, fmt = None, attr_mod = '', memsave = False,
                                    collapse_id = False, collapse_parent = True, fout_feature = None):
    ## functions to make list of id/parent to write
    if collapse_id: make_ids = lambda l: ';'.join(sorted(l))
    else: make_ids = lambda l: l
    if collapse_parent: make_parents = lambda l: ';'.join(sorted(l))
    else: make_parents = lambda l: l
    ## read into GFF object
    print("Reading files")
    if fmt is None:
        fmt = "BED" if gff_fname[-3:].upper() == "BED" else "GFF"
    gff = GFF(fname = gff_fname, fmt = fmt, attr_mod = attr_mod, memsave = memsave)
    ## get relationships + write
    print("Retrieving relationships")
    relationships = {} ## format: {<child>: {<set of parents>}}
    if fout_feature is not None:
        f_feature = open(fout_feature, "w+")
    with open(fout, "w+") as f:
        for entry in gff:
            ids = make_ids(entry.get_attr("ID", fmt = list))
            parents = make_parents(entry.get_attr("Parent", fmt = list))
            for id_str in ids:
                if id_str not in relationships:
                    relationships[id_str] = set()
                    if fout_feature is not None:
                        _ = f_feature.write(f"{id_str}\t{entry.feature}\n")
                for parent_str in parents:
                    ## skip if relationship has already been logged
                    if parent_str in relationships[id_str]:
                        continue
                    else:
                        relationships[id_str].update({parent_str})
                        _ = f.write(f"{id_str}\t{parent_str}\n")
    return

parser = argparse.ArgumentParser(description="extract unique parent-child relationships")
parser.add_argument("gff", type=str, help="path to GFF3/BED file")
parser.add_argument("fout", type=str, help="output file")
parser.add_argument("--memsave", help="memory-saving mode", action="store_true", default=False)
parser.add_argument("--fmt", type=str, help="input format [GFF/BED]")
parser.add_argument("--attr-mod", type=str, help="attribute modification")
parser.add_argument("--fout-feature", type=str, help="output file for feature types (optional)")
parser.add_argument("--collapse-id=", choices=("true", "false"), default="false",
                    dest="collapse_id", type=str, help="collapse ID attribute for same annotation")
parser.add_argument("--collapse-parent=", choices=("true", "false"), default="false",
                    dest="collapse_parent", type=str, help="collapse Parent attribute for same annotation")
args = parser.parse_args()

# print(args)

get_and_write_gff_relationships(
    args.gff, args.fout, fmt = args.fmt,
    attr_mod = args.attr_mod, memsave = args.memsave,
    collapse_id = (args.collapse_id == "true"),
    collapse_parent = (args.collapse_parent == "true"),
    fout_feature = args.fout_feature
)
