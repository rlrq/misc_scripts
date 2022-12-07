#!/usr/bin/python3
import argparse

import sys
# import tempfile
sys.path.append("/mnt/chaelab/rachelle/src")

from fasta_manip import dict_to_fasta

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description=("remove entries with repeated names"
                 " (only for files small enough to be read into memory)"))
parser.add_argument("fname", type=str, help="path to input FASTA file")
parser.add_argument("-i", help="modify file in-place", action="store_true", default=False)
parser.add_argument("-k", help="keep duplicates (append '-<number>' to name); not compatible with -r",
                    action="store_true", default=False)
parser.add_argument("-r", help=("remove duplicates ONLY IF sequence AND sequence ID are identical;"
                                " keep and append '-<number>' to name otherwise;"
                                " not compatible with -k"),
                    action="store_true", default=False)
parser.add_argument("--out", type=str, help="path to output FASTA file")
args = parser.parse_args()

if args.r and args.k:
    raise Exception("-k and -r are mutually exclusive")
elif args.out and args.i:
    raise Exception("--out and -i are mutually exclusive")
elif not args.out and not args.i:
    raise Exception("either --out or -i is required")

output = {}
count = {}
for record in SeqIO.parse(args.fname, "fasta"):
    if not record.id in output:
        count[record.id] = [record.id]
        output[record.id] = record.seq
    else:
        ## ignore record if seqs w/ duplicate ids are to be removed
        if not args.r and not args.k:
            continue
        ## ignore record if complete duplicate of existing record (same id + same seq)
        elif args.r:
            if any([output[seqid] == record.seq for seqid in count[record.id]]):
                continue
        ## store new seq
        new_seqid = f"{record.id}-{len(count[record.id])+1}"
        output[new_seqid] = record.seq
        count[record.id].append(new_seqid)

## write
fout = args.fname if args.i else args.out
dict_to_fasta(output, fout)
