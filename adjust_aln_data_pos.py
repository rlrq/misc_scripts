#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from Bio import SeqIO

def adjust_data_pos(f_data, f_aln, f_out, seqid, col_pos = 0, col_data = 1, offset_in = 0, offset_out = 0,
                    numeric_only = True, header = False, gap_char = '-'):
    seq = None
    for record in SeqIO.parse(f_aln, "fasta"):
        if record.description == seqid:
            seq = record.seq
            break
    ## exit if no sequences match the provided ID
    if seq is None:
        f_out.write('')
        print(f"Sequence ID '{seqid}' not found.")
        return
    ## skip header line
    if header:
        next(f_data)
    ## parse data
    data = {}
    for line in f_data:
        curr_data = line.replace('\n', '').split('\t')
        data[int(curr_data[col_pos])] = curr_data[col_data]
    ## adjust pos and write data
    def is_numeric(v):
        try: float(v)
        except: return False
        return True
    i_ungapped = offset_out
    for i, c in enumerate(seq):
        i_aln = i - offset_in
        if (c == gap_char
            or i_aln not in data
            or (numeric_only and not is_numeric(data[i_aln]))):
            continue
        i_ungapped += 1
        f_out.write('\t'.join(map(str, [i_ungapped, data[i_aln]])) + '\n')
    return

parser = argparse.ArgumentParser(description="prepare file for loading into PyMOL as b-factor; removes all positions not present in sequence with specified ID (--id) and converts positions in alignment back to ungapped positions in specified sequence")
parser.add_argument("f_data", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("--fout", nargs='?', type=argparse.FileType('w+'), default=sys.stdout)
parser.add_argument("--aln", type=argparse.FileType('r'), help="file of FASTA alignment")
parser.add_argument("--id", type=str, help="ID of sequence to filter data and adjust positions to")
parser.add_argument("--pos", type=int, help="index of column (1-index) containing position (1-index)", default=1)
parser.add_argument("--data", type=int, help="index of column (1-index) containing data values", default=2)
parser.add_argument("--keep-non-numeric",
                    help="keep all data values; if not raised, only numeric values are retained",
                    action="store_true", default=False)
parser.add_argument("--header",
                    help="raise flag if header is present; header will be removed",
                    action="store_true", default=False)
parser.add_argument("--gap-char", type=str, help="gap character (default='-')", default='-')
parser.add_argument("--offset-in", type=int,
                    help=("position offset in alignment;"
                          " (e.g. if '--offset-in 4', position 1 in data becomes position 5 in alignment)"),
                    default=0)
parser.add_argument("--offset-out", type=int,
                    help=("position offset in output file;"
                          " (e.g. if '--offset-out 4', position 1 (post-de-gapping) becomes position 5)"),
                    default=0)
args = parser.parse_args()

adjust_data_pos(args.f_data, args.aln, args.fout, args.id,
                col_pos = args.pos-1, col_data = args.data-1, ## convert 1-indexed column index to 0-index
                offset_in = args.offset_in, offset_out = args.offset_out,
                numeric_only = (not args.keep_non_numeric),
                header = args.header, gap_char = args.gap_char)
