#!/usr/bin/python3
import argparse

import regex as re
from Bio import SeqIO
from Bio import Seq
from Bio import Data


parser = argparse.ArgumentParser(description="find motifs in FASTA file")
parser.add_argument("fa_in", type=str, help="path to FASTA file")
parser.add_argument("fout", type=str, help="output file")
parser.add_argument("motif", type=str, help="comma-separated sequence motifs")
parser.add_argument("--allow-palindrome", action="store_true", default=False,
                    help="allow palindrome (i.e. search also for reverse of motif without complementation)")
parser.add_argument("--chrom-pattern", type=str, default="^[^:]+",
                    help=("regex pattern to get chromosome from sequence ID;"
                          " use if FASTA file was subsetted from original genome assembly file"))
parser.add_argument("--offset-pattern", type=str, default="(?<=:)\d+",
                    help=("regex pattern to get coordinate offset from sequence ID;"
                          " use if FASTA file was subsetted from original genome assembly file"))
parser.add_argument("--adjust-coords", help="include columns for adjusted coordinates", action="store_true", default=False)
parser.add_argument("--group-distance", type=int, default=10,
                    help=("distance (bp) within which to group motifs;"
                          " note that motifs that are adjacent but not overlapping are defined as 1bp apart"))
args = parser.parse_args()

fa_in = args.fa_in
fout = args.fout
motifs = args.motif.split(',')
allow_palindrome = args.allow_palindrome
adjust_coords = args.adjust_coords
chrom_pattern = args.chrom_pattern
offset_pattern = args.offset_pattern
group_threshold = args.group_distance

def expand_ambiguous_bases(seq):
    output = ''
    for base in seq:
        expanded_base = Data.IUPACData.ambiguous_dna_values.get(base.upper(), base)
        if len(expanded_base) > 1:
            expanded_base = '[' + expanded_base + ']'
        output += expanded_base
    return output

# fa_in = "/mnt/chaelab/rachelle/zzOtherzz/Zimin/mute/phytozome14/20251202/167/fasta/167.subset.fasta"
# fout = "/mnt/chaelab/rachelle/zzOtherzz/Zimin/mute/phytozome14/20251202/167/167.subset.motif_TGTCNN.tsv"
# motifs = ["TGTCNN"]
# group_threshold = 10 ## bp
# allow_palindrome = True

if allow_palindrome:
    motifs = motifs + [motif[-1::-1] for motif in motifs]

motifs_revcomp = [Seq.Seq(motif).reverse_complement() for motif in motifs]
motifs_regex_fwd = '|'.join([expand_ambiguous_bases(motif) for motif in motifs])
motifs_regex_rvs = '|'.join([expand_ambiguous_bases(motif) for motif in motifs_revcomp])

class MotifMatch():
    def __init__(self, match_object, strand):
        self.match = match_object
        self.strand = strand
        self.span = self.match.span()
        self.start = min(self.span)
        self.end = max(self.span)
        self.seq = self.match.group(0)
    def __repr__(self):
        return self.match.__repr__()

## get maches
matches = {record.id: ([MotifMatch(match, '+') for match in re.finditer(motifs_regex_fwd, str(record.seq))] +
                       [MotifMatch(match, '-') for match in re.finditer(motifs_regex_rvs, str(record.seq))])
           for record in SeqIO.parse(fa_in, "fasta")}

## group by proximity
all_groups = {}
for seqid, seq_matches in matches.items():
    if not seq_matches: continue
    sorted_matches = sorted(seq_matches, key = lambda m:m.match.span())
    groups = []
    current_group = [sorted_matches[0]]
    current_range = list(sorted_matches[0].span)
    for match in sorted_matches[1:]:
        if match.start - current_range[1] <= group_threshold:
            current_group.append(match)
            current_range = [min(current_range + [match.start]), max(current_range + [match.end])]
        else:
            groups.append(current_group)
            current_group = [match]
            current_range = list(match.span)
    groups.append(current_group)
    all_groups[seqid] = groups

## write (start & end are 1-indexed and inclusive)
header = ["seqid", "start", "end", "strand", "sequence", "group", "bp_from_last", "bp_to_next"]
if adjust_coords:
    header.extend(["true_chrom", "true_start", "true_end"]) ## true_start and true_end are based on OG coordinates
    get_chrom = lambda seqid:re.search(chrom_pattern, seqid).group(0)
    get_offset = lambda seqid:int(re.search(offset_pattern, seqid).group(0))
else:
    get_chrom = lambda seqid:seqid
    get_offset = lambda seqid:1

with open(fout, "w+") as f:
    _ = f.write('\t'.join(header) + '\n')
    for seqid, groups in all_groups.items():
        chrom = get_chrom(seqid)
        offset = get_offset(seqid)
        for i, group_matches in enumerate(groups):
            distances = {i: (group_matches[i+1].start - group_matches[i].end + 1) for i in range(len(group_matches)-1)}
            for j, match in enumerate(group_matches):
                to_write = [seqid, match.start+1, match.end, match.strand, match.seq,
                            i, distances.get(j-1, "NA"), distances.get(j, "NA")]
                if adjust_coords:
                    to_write.extend([chrom, match.start+offset, match.end+offset-1])
                _ = f.write('\t'.join(map(str, to_write)) + '\n')
