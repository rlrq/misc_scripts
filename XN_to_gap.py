#!/usr/bin/python3
import argparse

import sys
# sys.path.append("/mnt/chaelab/rachelle/src")

from Bio import SeqIO

def char_to_gap(record, char, gap_char = '-'):
    record.seq = record.seq.replace(char, gap_char)
    return record

def XN_to_gap(fin, fout, char, gap_char = '-'):
    records = SeqIO.parse(fin, "fasta")
    SeqIO.write((char_to_gap(record, char, gap_char = gap_char) for record in records),
                fout, "fasta")
    return

parser = argparse.ArgumentParser(description="replace character (like unknown amino acid X or nucleotide N) with gap in FASTA file")
parser.add_argument("fin", type=str, help="path to input FASTA file")
parser.add_argument("fout", type=str, help="path to output FASTA file")
parser.add_argument("char", type=str, help="character to replace with gap")
parser.add_argument("--gap-char", type=str, help="gap character", default='-')
args = parser.parse_args()

XN_to_gap(args.fin, args.fout, args.char, gap_char = args.gap_char)

# fin, fout, char = sys.argv[1:4]
# gap_char = '-' if len(sys.argv) < 5 else sys.argv[4]
# XN_to_gap(fin, fout, char, gap_char = gap_char)
