#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import fasta_manip

# def remove_stop_codon_and_write(fname, fout):
#     fasta_manip.remove_stop_codon(fname, fout)
#     return

parser = argparse.ArgumentParser(description="remove stop codons from FASTA file")
parser.add_argument("fin", type=str, help="path to input FASTA file")
parser.add_argument("fout", type=str, help="path to output FASTA file")
args = parser.parse_args()

fasta_manip.remove_stop_codon(args.fin, args.fout)
