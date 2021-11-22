#!/usr/bin/python3

import sys
sys.path.append("/mnt/chaelab/rachelle/src")
from fasta_manip import fasta_to_dict, dict_to_fasta

def translate_fasta(fin, fout, table = 1):
    seqs = fasta_to_dict(fin)
    seqs = {k: v.translate(table = table) for k, v in seqs.items()}
    dict_to_fasta(seqs, fout)
    return

fin, fout = sys.argv[1:3]
table = 1 if len(sys.argv) < 4 else sys.argv[3]

translate_fasta(fin, fout, table = table)
