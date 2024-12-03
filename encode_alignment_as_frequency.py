#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import basics
from Bio import SeqIO

## replace all characters (except gaps) in alignment with rank of frequency in that position
def encode_alignment_as_frequency(fin, fout, gap_char = '-', delim = ',', keep_description = False, index = 1):
    ## read sequences
    records = list(SeqIO.parse(fin, "fasta"))
    ## get longest sequence so we can iterate through all
    ## (ideally the sequences will be of the same length but you never know)
    max_length = max(len(r) for r in records)
    ## create list to store encoded output as matrix first
    output = [[] for i in range(len(records))]
    ## iterate through each position
    exclude_chars = ['-']
    for i in range(max_length):
        ## get characters at position
        chars = [r.seq[i] for r in records]
        ## encode as frequency rank
        chars_ranked = basics.order_to_dict(basics.order_elements_by_freq(chars, exclude = exclude_chars), index = index)
        ## store output in alignment matrix
        for j, c in enumerate(chars):
            output[j].append(chars_ranked.get(c, c))
    ## write encoded alignment
    make_record_id = (lambda record: record.description) if keep_description else (lambda record: record.id)
    with open(fout, "w+") as f:
        for i, encoded in enumerate(output):
            f.write(','.join(map(str, [make_record_id(records[i])] + encoded)) + '\n')
    return


parser = argparse.ArgumentParser(description="encode characters in FASTA file by position as rank of frequency at position")
parser.add_argument("fin", type=str, help="path to input FASTA file")
parser.add_argument("fout", type=str, help="path to output delimited file")
parser.add_argument("--keep-description", help="use full sequence description as ID in output", action="store_true", default=False)
parser.add_argument("--delim", type=str, help="delimiter of output file", default=',')
parser.add_argument("--gap-char", type=str, help="gap character", default='-')
parser.add_argument("--index", type=int, help="starting index for ranking", default=1)
args = parser.parse_args()

encode_alignment_as_frequency(args.fin, args.fout,
                              gap_char = args.gap_char, delim = args.delim,
                              keep_description = args.keep_description, index = args.index)

