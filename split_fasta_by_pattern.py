#!/usr/bin/python3
import argparse

import re
import os.path
from Bio import SeqIO

def split_fasta_by_pattern(fin, fout_function, split_function, append = False):
    ## split sequences into groups
    output = {}
    for seq_record in SeqIO.parse(fin, "fasta"):
        grp = split_function(seq_record)
        output[grp] = output.get(grp, [])
        output[grp].append(seq_record)
    ## write output file
    for grp, seq_records in output.items():
        fout = fout_function(grp)
        if append and os.path.exists(fout) and os.path.isfile(fout):
            existing_seqs = list(SeqIO.parse(fout, "fasta"))
        else:
            existing_seqs = []
        SeqIO.write(existing_seqs + seq_records, fout, "fasta")
    return

parser = argparse.ArgumentParser(description="Split FASTA file according to patterns of sequence names.\nTakes as input a FASTA file, a function to generate output file name, and a function that assigns groups to each FASTA entry.")
parser.add_argument("fin", type=str, help="path to input FASTA file for splitting")
parser.add_argument("fout_function", type=str,
                    help="function that accepts a str (group name) to generate output FASTA file path")
parser.add_argument("split_function", type=str,
                    help="function that accepts Bio.SeqRecord.SeqRecord object and outputs a group name (str)")
parser.add_argument("--append", action="store_true", default=False,
                    help="if output file already exists, append to it instead of overwriting")
# parser.add_argument("--acc", type=str, help="comma-separated profile accessions to subset",
#                     default = '')
# parser.add_argument("--acc-v", type=str,
#                     help="comma-separated profile accessions with version to subset",
#                     default = '')
# parser.add_argument("--exclude", type=str,
#                     help="comma-separated field names to exclude",
#                     default = '')
# parser.add_argument("--essential", action="store_true", default=False,
#                     help="write only essential fields")
args = parser.parse_args()

# ## example usage
# split_fasta_by_pattern("/path/to/file.fasta",
#                        lambda grp: f"/path/to/output/dir/prefix.{grp}.fasta",
#                        lambda seq_record: seq_record.id.split('|')[2],
#                        append = True)

split_fasta_by_pattern(args.fin, eval(args.fout_function), eval(args.split_function), append = args.append)
