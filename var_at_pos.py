#!/usr/bin/python3

"""
take FASTA alignment and outputs (1) tsv of counts at pos and/or (2) trimmed alignment
based on either (1) ungapped position in reference sequence or (2) position in alignment
"""

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from fasta_manip import adjusted_pos, fasta_to_dict
from basics import get_count_dict

## returns dict of {<ungapped ref pos>:<aln pos>}
def reference_pos_in_alignment(seq, *positions, index = 1):
    output = {}
    for pos in positions:
        ## adjusted_pos is 1-indexed
        ## output pos uses input index system
        output[pos] = adjusted_pos(seq, pos + (1-index)) - (1-index)
    return output

## extract characters in relevant position
## aln: dictionary of {<seqid>: <seq>}
def chars_in_pos(aln, pos, index = 1):
    return {k: seq[pos-index] for k, seq in aln.items()}

## extract characters in multiple positions
def chars_in_pos_multi(aln, *positions, index = 1):
    return {pos: chars_in_pos(aln, pos, index = index) for pos in positions}

## makes count dict of characters at each position
def chars_in_pos_summ(aln, *positions, index = 1):
    multi_pos_data = chars_in_pos_multi(aln, *positions, index = index)
    return {pos: get_count_dict([str(c) for c in data.values()]) for pos, data in multi_pos_data.items()}

## gets count dict of characters at each position AND of reference sequence
def get_chars_in_pos_summ(aln, ref_seqid, *positions, index = 1):
    all_seqs_summ = chars_in_pos_summ(aln, *positions, index = index)
    multi_pos_data = chars_in_pos_multi(aln, *positions, index = index)
    ref_seq_char = {pos: data[ref_seqid] for pos, data in multi_pos_data.items()}
    return all_seqs_summ, ref_seq_char

## main function
## written positions will be 1-indexed
def var_at_pos(fout, fasta, ref_seqid, *positions_ref_ungapped, index = 1):
    aln = fasta_to_dict(fasta)
    if ref_seqid not in aln:
        raise Exception(f"Reference sequence ID ('{ref_seqid}') not in alignment.")
    pos_ref_to_aln = reference_pos_in_alignment(aln[ref_seqid], *positions_ref_ungapped, index = index)
    all_seqs_summ, ref_seq_char = get_chars_in_pos_summ(aln, ref_seqid, *pos_ref_to_aln.values(), index = index)
    with open(fout, "w+") as f:
        _ = f.write('\t'.join(["pos.ref", "pos.aln", "char.ref", "char.counts"]) + '\n')
        for pos_ref, pos_aln in sorted(pos_ref_to_aln.items()):
            pos_ref_index1 = pos_ref + (1-index)
            pos_aln_index1 = pos_aln + (1-index)
            char_ref = ref_seq_char[pos_aln]
            char_summ = all_seqs_summ[pos_aln]
            sorted_char_summ = sorted(char_summ.items(), key = lambda e: (-e[1],e[0]))
            _ = f.write('\t'.join(map(str,
                                      [pos_ref_index1, pos_aln_index1, char_ref,
                                       ','.join(f"{char}:{count}" for char, count in sorted_char_summ)]
            )) + '\n')
    return

# ## PAP2 lin2023
# var_at_pos(
#     "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/variation/PAP2.Col-0_lin2023-pseudomolecule.pep.mafft.interfaceResidues.tsv",
#     "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/aln/PAP2.Col-0_lin2023-pseudomolecule.pep.mafft.fasta",
#     "Col-0_ref|AT3G61540|CDS|AT3G61540.1",
#     165,167,169,170,171,172,332,335,345,349,425,
#     index = 1
# )

# ## PAP2 1001genomes
# var_at_pos(
#     "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/variation/PAP2.1001genomes.pep.interfaceResidues.tsv",
#     "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/seq/PAP2.1001genomes.pep.fasta",
#     "6909|Col-0|AT3G61540|CDS|AT3G61540.1",
#     165,167,169,170,171,172,332,335,345,349,425,
#     index = 1
# )
