#!/usr/bin/python3
import argparse
import sys
import copy

sys.path.append("/mnt/chaelab/rachelle/src")

from index_fasta import IndexedFasta
from fasta_manip import dict_to_fasta

def remove_column(aln, gap_freq = 1, gap_char = '-'):
    """
    Arguments:
        aln (dict/IndexedFasta): alignment
        gap_freq (float): maximum gap frequency to retain a column (default=1);
            columns with greater than 'gap_freq' frequency of gaps
            will be removed from the alignment;
            should be between 0 and 1.
        gap_char (str): gap character (default='-')
    
    Returns:
        dict: of trimmed sequences
    """
    ## get some constants
    seq_ids = tuple(aln.keys())
    num_seq = len(seq_ids)
    aln_len = len(aln[seq_ids[0]])
    empty_seq = aln[seq_ids[0]][:0]
    progress_increment = set(int(i*aln_len/10) for i in range(10))
    ## start trimming
    output = {}
    for i in range(aln_len):
        if i in progress_increment:
            print(f"Processing position {i} of {num_seq}")
        col = tuple(aln[seq_id][i] for seq_id in seq_ids)
        if (len([c for c in col if c == gap_char])/num_seq <= gap_freq):
            for i, seq_id in enumerate(seq_ids):
                output[seq_id] = output.get(seq_id, empty_seq) + col[i]
    return output


parser = argparse.ArgumentParser(description="trim FASTA file of alignment by gap frequency per column")
parser.add_argument("--in", type=str, help="path to input FASTA alignment file",
                    dest = "fin")
parser.add_argument("--out", type=str, help="path to output (trimmed) FASTA alignment file",
                    dest = "fout")
parser.add_argument("--gap-char", type=str, help="gap character",
                    default = '-')
parser.add_argument("--gap-freq", type=float, help="maximum gap frequency for retained columns",
                    default = 1)
args = parser.parse_args()

aln = IndexedFasta(args.fin)
aln_trimmed = remove_column(aln, gap_freq = args.gap_freq, gap_char = args.gap_char)
dict_to_fasta(aln_trimmed, args.fout)
