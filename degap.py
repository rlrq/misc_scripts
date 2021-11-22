#!/usr/bin/python3

import sys
sys.path.append("/mnt/chaelab/rachelle/src")
from fasta_manip import fasta_to_dict, dict_to_fasta

def degap_fasta(fin, fout, gap_char = '-'):
    seqs = fasta_to_dict(fin)
    seqs = {k: str(v).replace(gap_char, '') for k, v in seqs.items()}
    dict_to_fasta(seqs, fout)
    return

fin, fout = sys.argv[1:3]
gap_char = '-' if len(sys.argv) < 4 else sys.argv[3]

degap_fasta(fin, fout, gap_char = gap_char)
