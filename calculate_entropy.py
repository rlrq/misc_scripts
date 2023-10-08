#!/usr/bin/python3
import argparse
from scipy.stats import entropy

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from basics import get_count_dict
from fasta_manip import fasta_to_dict

## get discrete distribution of raw sampled data (e.g. ['a', 'a', 'a', 'b', 'b', 'c', 'c'] -> [3, 2, 2])
## not order preserving
def make_discrete_distr(array):
    counts = tuple(get_count_dict(array).values())
    total = sum(counts)
    return tuple(count/total for count in counts)

## calculate sitewise entropy for sequences, returns list of entropy values same length as alignment
def entropy_sitewise(seqs, ignore_gaps = True, gap_char = '-', stop_codon = '*'):
    aln_len = len(seqs[0])
    output = [None] * aln_len
    for i in range(aln_len):
        aa = [seq[i] for seq in seqs if seq[i] != stop_codon]
        if ignore_gaps:
            aa = [a for a in aa if a != gap_char]
        output[i] = entropy(make_discrete_distr(aa))
    return output

## calculate entropy using all sequences in alignment
def entropy_from_all(f_aln, fout, **kwargs):
    seqs = tuple(fasta_to_dict(f_aln).values())
    ent = entropy_sitewise(seqs, **kwargs)
    with open(fout, 'w+') as f:
        f.write('\t'.join(["grp", "pos", "H"]) + '\n')
        for i, val in enumerate(ent):
            f.write('\t'.join(["master", str(i+1), str(val)]) + '\n')
    return

## read grouping file (format: <grp ID>\t<seq ID>) (no header)
def read_groups(fname):
    output = {}
    with open(fname, 'r') as f:
        for line in f.readlines():
            grp, seqid = line.replace('\n', '').split('\t')
            output[grp] = output.get(grp, set()).union({seqid})
    return output

def entropy_by_group(f_aln, fout, f_grp, **kwargs):
    ## parse seqs
    seqs = fasta_to_dict(f_aln)
    ## parse groups
    grps = read_groups(f_grp)
    ## open output file
    with open(fout, 'w+') as f:
        f.write('\t'.join(["grp", "pos", "H"]) + '\n')
        for grpid, seqids in grps.items():
            ## calculate entropy
            ent = entropy_sitewise([seq for seqid, seq in seqs.items() if seqid in seqids], **kwargs)
            ## write
            for i, val in enumerate(ent):
                f.write('\t'.join([grpid, str(i+1), str(val)]) + '\n')
    return


parser = argparse.ArgumentParser(description="calculate entropy from alignment")
parser.add_argument("f_aln", type=str, help="path to FASTA alignment file")
parser.add_argument("fout", type=str, help="output file")
parser.add_argument("f_grp", type=str, help="taxon grouping file; if provided, entropy will be calculated only within taxon", nargs='?', default=None)
parser.add_argument("--stop-codon", type=str, help="stop codon character (default='*')", default='*')
parser.add_argument("--gap", type=str, help="gap character (default='-')", default='-')
parser.add_argument("--include-gaps", help="include gap characters in entroy calculation", action="store_true", default=False)
args = parser.parse_args()

if args.f_grp is None:
    entropy_from_all(args.f_aln, args.fout,
                     gap_char = args.gap,
                     stop_codon = args.stop_codon,
                     ignore_gaps = (not args.include_gaps))
else:
    entropy_by_group(args.f_aln, args.fout, args.f_grp,
                     gap_char = args.gap,
                     stop_codon = args.stop_codon,
                     ignore_gaps = (not args.include_gaps))


