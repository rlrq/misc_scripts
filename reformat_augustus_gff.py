#!/usr/bin/python3

import argparse

## reformat augustus output (adjust positions and remove comment lines)
def reformat_augustus(gff_in, gff_out, prefix = None, remove_comments = True, subsetted_sequence = False):
    with open(gff_in, 'r') as fin:
        with open(gff_out, 'w+') as fout:
            if remove_comments:
                fout.write("##gff-version 3\n")
            for line in fin:
                if line[0] == '#':
                    if not remove_comments:
                        fout.write(line)
                    continue
                line = line.split('\t')
                ## determine chrom, start, and end
                if subsetted_sequence:
                    chrom, coords = line[0].split(':')
                    offset = int(coords.split('..')[0])
                    start = int(line[3])+offset-1
                    end = int(line[4])+offset-1
                else:
                    chrom = line[0]
                    start = line[3]
                    end = line[4]
                if prefix is not None:
                    attributes = line[-1].replace("ID=", f"ID={prefix}").replace("Parent=", f"Parent={prefix}")
                else:
                    attributes = line[-1]
                fout.write('\t'.join(map(str,
                                         [chrom] + line[1:3] +
                                         [start, end] +
                                         line[5:-1] +
                                         [attributes])))
    return

parser = argparse.ArgumentParser(description="reformat gff file output by Augustus")
parser.add_argument("--in", type=str, help="path to input gff file for reformatting",
                    dest = "fin")
parser.add_argument("--out", type=str, help="path to output (reformatted) gff file",
                    dest = "fout")
parser.add_argument("--keep-comments", type=str, help="retain all comment lines",
                    dest = "keep_comments")
parser.add_argument("--prefix", type=str, help="feature ID prefix",
                    dest = "prefix", default = None)
parser.add_argument("--subsetted-sequence",
                    help=("raise this flag if sequence are subsetted and coords are to be adjusted back to original,"
                          " where the format in the chromosome field is '<original seqid>:<subset start>..<subset end>',"
                          " where subset coordinates are 1-index and end-inclusive"),
                    dest = "subsetted_sequence", action = "store_true", default = False)
args = parser.parse_args()

reformat_augustus(args.fin, args.fout, prefix = args.prefix,
                  remove_comments = (not args.keep_comments),
                  subsetted_sequence = args.subsetted_sequence)
