#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from fasta_manip import adjusted_pos, fasta_to_dict
from range_manip import Range, Ranges

class Feature:
    def __init__(self, seqid, feature, start, end):
        self.seqid = seqid
        self.feature = feature
        self.range = Range(start, end)
        self.adjusted_range = {"complete": None,
                               "disjoint": None}
        self.start = start
        self.end = end
        return
    ## returns Ranges of Range (list of list)
    def range_in_aln(self, gapped_seq, complete = False, use_cache = False):
        range_type = "complete" if complete else "disjoint"
        if not use_cache or self.adjusted_range[range_type] is None:
            self.adjusted_range[range_type] = gapped_seq.adjusted_range(self.range,
                                                                        subtract_gaps = (not complete))
        return self.adjusted_range[range_type]

class GappedSeq:
    def __init__(self, seqid, seq):
        self.seqid = seqid
        self.seq = seq
        self.ranges = self._range()
        return
    ## parse self.seq to get range(s) of non-gap pos
    ## returns Ranges object
    def _range(self):
        ranges = Ranges()
        curr_start = None
        ## get ranges
        for i, c in enumerate(self.seq):
            ## start range
            if curr_start is None and c != '-':
                curr_start = i
            ## terminate range
            elif curr_start is not None and c == '-':
                ranges.append(Range(curr_start, i))
                curr_start = None
        ## check for unterminated range at end of loop
        if curr_start is not None:
            ranges.append(Range(curr_start, i+1))
        return ranges
    ## return gapped range(s) (list of tuple ranges) from ungapped range (start is inclusive, end is exclusive, both 0-index)
    def adjusted_range(self, range_obj, subtract_gaps = True):
        gapped_start = adjusted_pos(self.seq, range_obj.start + 1) - 1
        gapped_end = adjusted_pos(self.seq, range_obj.end)
        if not subtract_gaps:
            return Ranges([Range(gapped_start, gapped_end)])
        else:
            return GappedSeq("tmp", self.seq[gapped_start:gapped_end]).ranges.offset(gapped_start)
        
## file should have the following columns (with header row): seqid, feature_name, start, end
## returns list of Feature objects
## index and end_exclusive refer to attributes of ranges in fname
def parse_ungapped_ranges(fname, index = 1, end_exclusive = False):
    data = []
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            tmp = line.rstrip('\n').split('\t')
            tmp[2] = int(tmp[2]) - index
            tmp[3] = int(tmp[3]) - index + (0 if end_exclusive else 1)
            data.append(Feature(*tmp))
    return data

def parse_ungapped_ranges_from_aln(fname):
    data = []
    from Bio import SeqIO
    for record in SeqIO.parse(fname, "fasta"):
        data.append(Feature(record.description, "master", 0, len(record.seq.replace('-', ''))))
    return data

## convert to range in alignment, returns list of [(Feature, [(start, end,), (start, end,)])]
## features: list of [Feature]
## aln: dict of {seqid: <gapped seq>}
def adjust_range_per_seq(features, aln, complete = False, use_cache = False):
    output = []
    for feature in features:
        if feature.seqid not in aln:
            output.append((feature, Ranges()))
        else:
            output.append(
                (feature,
                 feature.range_in_aln(aln[feature.seqid],
                                      complete = complete,
                                      use_cache = use_cache)
                )
            )
    return output


## combine overlapping/adjacent ranges of same feature
## features: list of [Feature]
## aln: dict of {seqid: <gapped seq>}
def adjust_range_per_feature(features, aln, use_cache = False):
    feature_ranges = {}
    for feature in features:
        if feature.seqid not in aln:
            continue
        else:
            feature_ranges[feature.feature] = feature_ranges.get(feature.feature, Ranges).union(
                feature.range_in_aln(aln[feature.seqid], complete = False, use_cache = use_cache)
            )
    return feature_ranges

def adjust_ranges(fname, aln_fname, fout_by_seq = None, fout_by_feature = None, complete = False,
                  index = 1, end_exclusive = False):
    aln = {seqid: GappedSeq(seqid, seq) for seqid, seq in fasta_to_dict(aln_fname).items()}
    if fname is not None:
        features = parse_ungapped_ranges(fname)
    else:
        features = parse_ungapped_ranges_from_aln(aln_fname)
    if fout_by_seq:
        data = adjust_range_per_seq(features, aln, complete = complete)
        data.sort(key = lambda entry: (entry[0].seqid, entry[0].feature, entry[0].start, entry[0].end))
        with open(fout_by_seq, 'w+') as f:
            f.write('\t'.join(["seqid", "feature", "ranges"]) + '\n')
            for entry in data:
                feature, ranges = entry
                f.write('\t'.join([feature.seqid,
                                   feature.feature,
                                   ranges.to_string(index = index, end_exclusive = end_exclusive)]) +
                        '\n')
    if fout_by_feature:
        data = adjust_range_per_feature(features, aln, use_cache = True)
        with open(fout_by_feature, 'w+') as f:
            f.write('\t'.join(["feature", "ranges"]) + '\n')
            for feature in sorted(data.keys()):
                f.write('\t'.join([feature, data[feature].to_string(index = index, end_exclusive = end_exclusive)]) + '\n')
    return


parser = argparse.ArgumentParser(description="get ranges in alignment")
parser.add_argument("aln", type=str, help="path to FASTA alignment file")
parser.add_argument("f_ranges", type=str, help="path to file containing ungapped ranges (<seqid>\t<range; 1-index, end-inclusive>)", nargs='?', default=None)
parser.add_argument("--fout-seq", type=str, help="output file of adjusted ranges by sequence", default=None)
parser.add_argument("--fout-feature", type=str, help="output file of adjusted ranges by feature", default=None)
parser.add_argument("--end-exclusive", help="write end-exclusive ranges", action="store_true", default=False)
parser.add_argument("--index", type=int, help="write ranges with index starting at <val>", default=1)
parser.add_argument("--complete", help="write ranges inclusive of intervening gaps", action="store_true", default=False)
args = parser.parse_args()

adjust_ranges(args.f_ranges, args.aln,
              fout_by_seq = args.fout_seq,
              fout_by_feature = args.fout_feature,
              complete = args.complete,
              index = args.index,
              end_exclusive = args.end_exclusive)
