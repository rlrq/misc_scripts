#!/usr/bin/python3

import argparse
import ast
import numpy as np

parser = argparse.ArgumentParser(description = "summarise computeMatrix output row-wise")
parser.add_argument("input", metavar="input.matrix", help="Input matrix file")
parser.add_argument("output", metavar="output.txt", help="Output txt file")
parser.add_argument("--gz", help="input is gzipped", default=False, action="store_true")
parser.add_argument("--stat", help="summary statistic", default="mean", type=str,
                    choices=["mean", "median", "min", "max", "sum", "std"])
parser.add_argument("--multi-sample", help="calculate metrics separately for each sample",
                    default=False, action="store_true")
parser.add_argument("--nan-as-zero", help="treat NaN as 0", default=False, action="store_true")
args = parser.parse_args()

def make_summarise(stat, nan_as_zero = False):
    ## pre-process any NaNs
    if nan_as_zero: process = lambda vals: [(0 if np.isnan(x) else x) for x in vals]
    else: process = (lambda vals: vals)
    ## get metric function
    if stat == "mean": calc = lambda vals: np.nanmean(vals)
    elif stat == "median": calc = lambda vals: np.nanmedian(vals)
    elif stat == "min": calc = lambda vals: min(vals)
    elif stat == "max": calc = lambda vals: max(vals)
    elif stat == "sum": calc = lambda vals: np.nansum(vals)
    elif stat == "std": calc = lambda vals: np.nanstd(vals)
    return lambda vals: calc(process(vals))

## make function to calculate stats
summarise = make_summarise(args.stat, args.nan_as_zero)

def process_file(fin_handler):
    start = fin.tell()
    l1 = fin.readline()
    ## get expected # of columns per sample by reading header, if header exists
    if l1[0] == '@':
        header = ast.literal_eval(l1[1:].rstrip().replace(":false,",":False,").replace(":null,",":None,").replace(":true,",":True,"))
        sample_boundaries = header["sample_boundaries"]
        num_samples = len(header["sample_labels"])
        if not args.multi_sample:
            sample_boundaries = [0, sample_boundaries[-1]]
            num_samples = 1
    else:
        sample_boundaries = [0, len(l1.split('\t')) - 6] ## first 6 columns identify the feature
        num_samples = 1
        fin.seek(start) ## reset position to start of file so we don't skip processing this first line
    ## open output file
    with open(args.output, 'w+') as fout:
        ## write header
        header["sample_boundaries"] = list(range(num_samples+1))
        _ = fout.write('@' + str(header) + '\n')
        ## process
        for line in fin:
            line_split = line.rstrip().split('\t')
            ## first 6 columns identify the feature
            feature_meta = line_split[:6]
            ## next however many columns contain values to process
            values = [float(x) for x in line_split[6:]]
            ## split values by sample according to sample_boundaries
            values_by_sample = [values[sample_boundaries[i]:sample_boundaries[i+1]] for i in range(num_samples)]
            ## calculate stats
            stats_by_sample = [summarise(vals) for vals in values_by_sample]
            ## write
            _ = fout.write('\t'.join(feature_meta + [str(x) for x in stats_by_sample]) + '\n')
    ## finish
    return

## process
if args.gz or args.input[-3:] == ".gz":
    import gzip
    with gzip.open(args.input, 'rt') as fin:
        process_file(fin)
else:
    with open(args.input, 'r') as fin:
        process_file(fin)
## finish
