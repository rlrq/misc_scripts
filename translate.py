#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")
from fasta_manip import fasta_to_dict, dict_to_fasta

def translate_fasta(fin, fout, table = 1, no_stop = False):
    seqs = fasta_to_dict(fin)
    try:
        table = int(table)
    except ValueError:
        table = table
    seqs = {k: v.translate(table = table) for k, v in seqs.items()}
    if no_stop:
        seqs = {k: v.replace('*', '') for k, v in seqs.items()}
    dict_to_fasta(seqs, fout)
    return

# fin, fout = sys.argv[1:3]
# table = 1 if len(sys.argv) < 4 else sys.argv[3]

parser = argparse.ArgumentParser(description="translate nucleotide FASTA file")
parser.add_argument("fin", type=str, help="path to input FASTA file")
parser.add_argument("fout", type=str, help="path to output FASTA file")
parser.add_argument("table", type=str, help="codon table (per Biopython)", nargs='?', default=1)
parser.add_argument("--no-stop", help="exclude stop codon(s)", action="store_true", default=False)
args = parser.parse_args()

translate_fasta(args.fin, args.fout, table = args.table, no_stop = args.no_stop)
