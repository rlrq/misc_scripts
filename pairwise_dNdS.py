#!/usr/bin/python3
import argparse
import re
import dnds

from tqdm import tqdm

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from fasta_manip import fasta_to_dict

def keep_aligned_pos(s1, s2):
    s1_out = s1[0:0]
    s2_out = s2[0:0]
    for i in range(len(s1)):
        if s1[i] != '-' and s2[i] != '-':
            s1_out += s1[i]
            s2_out += s2[i]
    return s1_out, s2_out

# def dn_ds(s1, s2):
#     ## if no sequence or sequences are identical
#     if not s1 or s1 == s2:
#         return "NA"
#     try:
#         return dnds.dnds(str(s1).upper(), str(s2).upper())
#     except:
#         s1_pep = s1.translate()
#         s2_pep = s2.translate()
#         if all([s1_pep[i] != s2_pep[i]

## window = number of codons
def sliding_dn_ds(s1, s2, window = 10):
    window_nt = window * 3
    output = []
    for i in range(window_nt, len(s1) + 3, 3):
        s1_part, s2_part = keep_aligned_pos(s1[i-window_nt:i], s2[i-window_nt:i])
        try:
            output.append(dnds.dnds(str(s1_part).upper(), str(s2_part).upper()))
        except:
            output.append("NA")
    return ','.join(map(str, output))

def constant_dn_ds(s1, s2):
    s1, s2 = keep_aligned_pos(s1, s2)
    return dnds.dnds(str(s1).upper(), str(s2).upper())

def apply_pairwise(fasta, fout, function, metric_name = "value", **kwargs):
    seqs = fasta_to_dict(fasta)
    seq_ids = tuple(seqs.keys())
    with open(fout, "w+") as f:
        f.write(f"seq_1\tseq_2\t{metric_name}\n")
        for i in tqdm(range(len(seq_ids))):
            for j in range(len(seq_ids)):
                s1_id = seq_ids[i]
                s2_id = seq_ids[j]
                try:
                    val = function(seqs[s1_id], seqs[s2_id], **kwargs)
                except:
                    val = "NA"
                f.write(f"{s1_id}\t{s2_id}\t{val}\n")
    return

parser = argparse.ArgumentParser(description="subset GFF file")
parser.add_argument("fa_in", type=str, help="path to FASTA codon alignment")
parser.add_argument("fout", type=str, help="path to output TSV file")
parser.add_argument("--window", type=int, help="window size (number of codons); if not provided, constant dN/dS will be calculated", default=0)
args = parser.parse_args()

apply_pairwise(args.fa_in,
               args.fout,
               (constant_dn_ds if args.window == 0 else sliding_dn_ds),
               window = args.window,
               metric_name = "dnds")

# fa_in, fout = sys.argv[1:3]
# apply_pairwise(fa_in, fout, constant_dn_ds, metric_name = "dnds")

