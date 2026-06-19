#!/home/rachelle/miniconda3/envs/py3.13/bin/python
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import os
import itertools
import pandas as pd
from Bio import SeqIO
from plotnine import ggplot, aes
import plotnine as ggp

## note: for each x-axis sequence, a new directory will be created within dirout, named according to sequence ID
## for each y-axis sequence, a tsv file of score matrix + plot will be created within the x-axis' folder, named also according to sequence ID

parser = argparse.ArgumentParser(description="generate dot plot based on identity")
parser.add_argument("x", type=str, help="path to FASTA file for sequence(s) on x axis")
parser.add_argument("y", type=str, help="path to FASTA file for sequence(s) on y axis")
parser.add_argument("dirout", type=str, help="output directory")
parser.add_argument("--window", type=int, help="window size (bp)", default=4)
parser.add_argument("--offset", type=int, help="offset between windows", default=2)
parser.add_argument("--matrix-only", help="skip plot generation", action="store_true", default=False)
parser.add_argument("--min-pid", type=float, help="minimum % identity for plotting (default=50)")
parser.add_argument("--max-pid", type=float, help="minimum % identity for maximum plot intensity (default=100)")
parser.add_argument("--min-id", type=int, help=("minimum # identity in window for plotting"
                                                " (default=half window size)"))
parser.add_argument("--max-id", type=int, help=("minimum # identity in window for maximum plot intensity"
                                                " (default=window size)"))
parser.add_argument("--min-colour", type=str, help="colour for minimum plotted % identity", default="black")
parser.add_argument("--max-colour", type=str, help="colour for maximum plotted % identity", default="black")
parser.add_argument("--exclude-revcomp", help="only score matched strands",
                    action="store_true", default=False)
parser.add_argument("--coord-free", help="allow axes to scale independently",
                    action="store_true", default=False)
parser.add_argument("--binary", help="binarise colours based on whether identity is less than min_id",
                    action="store_true", default=False)
args = parser.parse_args()

# class Tmp():
#     def __init__(self):
#         return

# args = Tmp()
# args.window = 4
# args.offset = 2
# args.min_id = 50
# args.max_id = 100
# args.min_colour = "black"
# args.max_colour = "black"
# args.include_revcomp = True
# args.exclude_revcomp = False

# recordX = SeqIO.read("/mnt/chaelab/rachelle/zzOtherzz/Zimin/pin3/data/fasta/pin3.TAIR10.intron3.fasta", "fasta")
# recordX = SeqIO.read("/mnt/chaelab/rachelle/zzOtherzz/Zimin/pin3/data/fasta/flanked_AT1G70940.fasta", "fasta")
# recordY = SeqIO.read("/mnt/chaelab/rachelle/zzOtherzz/Zimin/pin3/data/fasta/pin4.TAIR10.intron3.fasta", "fasta")
# recordY = SeqIO.read("/mnt/chaelab/rachelle/zzOtherzz/Zimin/pin3/data/fasta/pin7.TAIR10.intron3.fasta", "fasta")

# class DotSeq():
#     def __init__(self, seq_record, window_size = args.window, offset = args.offset):
#         self.record = seq_record
#         self.seq = self.record.seq
#         self.window_seqs = []
#         self.window_hash = []
#         self.score_hash = {}
#         self.hash_windows(window_size, offset)
#         return
#     ## create self.windows and self.window_hash
#     def hash_windows(self, window_size, offset):
#         tmp_window_seqs = {}
#         for i in range(0, len(self.seq), offset):
#             ## get sequence within window
#             window_seq = self.seq[i:i+window_size]
#             ## add sequence to dict it not encountered before
#             if window_seq not in tmp_window_seqs:
#                 tmp_window_seqs[window_seq] = len(self.windows)
#             ## store index of sequence (index in self.window_hash is index of window)
#             self.window_hash.append(tmp_window_seqs[window_seq])
#         ## convert tmp_window_seqs to list
#         self.window_seqs = [0 for i in range(len(tmp_window_seqs))]
#         for seq, index in tmp_window_seqs.items():
#             self.window_seqs[index] = seq
#         return
#     ## create (temporary) dict for scores of windows between DotSeq objects
#     def hash_scores(self, other):
#         for window_seq in self.window_seqs
#         for other_i, other_seq in other.window_seqs

def id_count(seqA, seqB, include_revcomp = False):
    match_strands = sum([seqA[i] == seqB[i] for i in range(min(len(seqA), len(seqB)))])
    if include_revcomp:
        rev_strands = sum([seqA[i] == seqB[len(seqB)-i-1] for i in range(min(len(seqA), len(seqB)))])
        return max(match_strands, rev_strands)
    return match_strands

class DotMatrix():
    def __init__(self, recordX, recordY,
                 window = args.window, offset = args.offset, include_revcomp = (not args.exclude_revcomp)):
        self.recX = recordX
        self.recY = recordY
        self.matrix = None
        self.plot = None
        self.window = window
        self.offset = offset
        self.include_revcomp = include_revcomp
        self.make_matrix()
    def make_matrix(self, window = None, offset = None, include_revcomp = None):
        ## set parameters
        window = self.window if window is None else window
        offset = self.offset if offset is None else offset
        include_revcomp = self.include_revcomp if include_revcomp is None else include_revcomp
        self.window = window
        self.offset = offset
        self.include_revcomp = include_revcomp
        seqX = self.recX.seq
        seqY = self.recY.seq
        ## start scoring identity
        windows_Y = {}
        output = []
        for Yi in range(0, len(seqY), offset):
            winY = seqY[Yi:Yi+window]
            ## if window sequence is identical to a previous window, duplicate the scores for that window and skip to next window
            if winY in windows_Y:
                output.append(output[windows_Y[winY]])
                continue
            windows_Y[winY] = len(output) ## hash the new window sequence
            scoreY = [] ## create list for scores for this window
            ## score each window in X against Y
            windows_X = {}
            for Xi in range(0, len(seqX), offset):
                winX = seqX[Xi:Xi+window]
                if winX in windows_X:
                    scoreY.append(scoreY[windows_X[winX]])
                    continue
                scoreY.append(id_count(winY, winX, include_revcomp = include_revcomp))
            output.append(scoreY)
        self.matrix = output
        return
    def make_prefix(self):
        return f"{self.recY.id}.id.window-{self.window}_offset-{self.offset}_includeRevcomp-{self.include_revcomp}"
    def write_matrix(self, dirout, sep = ','):
        dirout_X = os.path.join(dirout, self.recX.id)
        os.makedirs(dirout_X, exist_ok = True)
        with open(os.path.join(dirout_X, self.make_prefix() + ".mat"), "w+") as f:
            for row in self.matrix:
                _ = f.write(sep.join(map(str,row)) + '\n')
        return
    def long_dataframe(self):
        df_long = pd.DataFrame(
            data = {"Xbin": (list(range(len(self.matrix[0])))*len(self.matrix)),
                    "Ybin": itertools.chain(*[[i]*len(self.matrix[0]) for i in range(len(self.matrix))]),
                    "val": itertools.chain(*self.matrix)})
        return df_long


# # mat_long = [list(range(len(dot_matrix.matrix[0])))*]
# dot_matrix = DotMatrix(recordX, recordY)
# df_long = dot_matrix.long_dataframe()


## calculate min and max
max_id = args.max_id if args.max_id is not None else (
    args.window if args.max_pid is None else round(args.window*(args.max_id/100)))
min_id = args.min_id if args.min_id is not None else (
    round(args.window/2) if args.min_pid is None else round(args.window*(args.min_id/100)))

def plot_dotplot(dot_matrix, min_id, max_id, coord_fixed = (not args.coord_free)):
    p = ggplot(dot_matrix.long_dataframe()) + \
        aes(x = "Xbin", y = "Ybin", fill = "val") + \
        ggp.geom_tile() + ggp.ylab(dot_matrix.recY.id) + ggp.xlab(dot_matrix.recX.id) + \
        ggp.scale_fill_cmap(limits = [min_id, max_id],
                            breaks = list(range(min_id, max_id+1))) + \
        ggp.scale_y_continuous(expand = [0,0]) + \
        ggp.scale_x_continuous(expand = [0,0])
    if coord_fixed:
        p = p + ggp.coord_fixed()
    # print(p) ## print (and p.show()) doesn't work since we don't have elpy
    return p

def plot_dotplot_binary(dot_matrix, min_id, coord_fixed = (not args.coord_free)):
    df = dot_matrix.long_dataframe()
    df = df.loc[(df["val"]>=min_id)]
    p = ggplot(df) + \
        aes(x = "Xbin", y = "Ybin") + \
        ggp.geom_tile(fill = args.max_colour) + \
        ggp.ylab(dot_matrix.recY.id) + \
        ggp.xlab(dot_matrix.recX.id) + \
        ggp.scale_y_continuous(expand = [0,0]) + \
        ggp.scale_x_continuous(expand = [0,0]) + \
        ggp.theme(panel_background = ggp.element_rect(fill = "white", colour = "black"))
    return p

for recordX in SeqIO.parse(args.x, "fasta"):
    print("X:", recordX.id)
    for recordY in SeqIO.parse(args.y, "fasta"):
        print("Y:", recordY.id)
        dot_matrix = DotMatrix(recordX, recordY)
        dot_matrix.write_matrix(os.path.realpath(args.dirout))
        if args.binary:
            p = plot_dotplot_binary(dot_matrix, min_id = min_id)
            p.save(os.path.join(os.path.realpath(args.dirout), recordX.id, dot_matrix.make_prefix() + ".binary.png"))
        else:
            print(min_id)
            p = plot_dotplot(dot_matrix, min_id = min_id, max_id = max_id)
            p.save(os.path.join(os.path.realpath(args.dirout), recordX.id, dot_matrix.make_prefix() + ".png"))

# def make_id_count_matrix(seqX, seqY, window, offset):
#     output = []
#     windows_Y = {}
#     for Yi in range(0, len(seqY), offset):
#         winY = seqY[Yi:Yi+window]
#         ## if window sequence is identical to a previous window, duplicate the scores for that window and skip to next window
#         if winY in windows_Y:
#             output.append(output[windows_Y[winY]])
#             continue
#         windows_Y[winY] = len(output)+1 ## hash the new window sequence
#         scoreY = [] ## create list for scores for this window
#         ## score each window in X against Y
#         windows_X = {}
#         for Xi in range(0, len(seqX), offset):
#             winX = seqX[Xi:Xi+window]
#             if winX in windows_X:
#                 scoreY.append(scoreY[windows_X[winX]])
#                 continue
#             scoreY.append(id_count(winY, winX))
#         output.append(scoreY)
#     return output

# def write_matrix(matrix, recordX, recordY, dirout, sep = '\t'):
#     dirout_X = os.path.join(dirout, recordX.id)
#     os.makedirs(dirout_X, exist_ok = True)
#     with open(os.path.join(dirout_X, recordY.id+".id.mat"), "w+"):
#         for row in matrix:
#             _ = f.write(sep.join(map(str,row)) + '\n')
#     return

# def plot_dotplot(matrix, recordX, recordY):
#     id_mat = id_count_matrix(recordX.seq, recordY.seq)
    
