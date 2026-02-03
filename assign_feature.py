#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

priority_order = [
    "CDS", "three_prime_UTR", "five_prime_UTR", "exon",
    "mRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "miRNA", "ncRNA", "gene",
    "pseudogenic_exon", "pseudogenic_transcript", "pseudogene",
    "mRNA_TE_gene", "transposable_element_gene"
]
priority_last = []
ignore = [ "protein", "chromosome" ]

parser = argparse.ArgumentParser(description="assign feature to range in BED file")
parser.add_argument("query", type=str, help="path to query BED file")
parser.add_argument("fout", type=str, help="path to output file (fmt: <query BED columns>\t<feature>")
parser.add_argument("-b", "--bed", type=str, nargs='+', help="path to database BED file of features")
parser.add_argument("-g", "--gff", type=str, nargs='+', help="path to database GFF3 file of features")
parser.add_argument("--ob", "--out-bed", type=str, help="path of file to write intersection w/ BED files")
parser.add_argument("--og", "--out-gff", type=str, help="path of file to write intersection w/ GFF3 files")
parser.add_argument("--oi", "--out-intersect", type=str, help="path of file to write intersection w/ all files")
parser.add_argument("--priority", type=str,
                    help=( "comma-separated features in order of most-to-least priority to prepend to priority list."
                           " Takes priority over --priority-last."))
parser.add_argument("--priority-last" type=str,
                    help="comma-separated features in order of most-to-least priority to append to priority list")
parser.add_argument("--report-only", type=str,
                    help=( "comma-separated features to report; all others will be collased as 'others'."
                           " Takes priority over --ignore." ))
parser.add_argument("--ignore", type=str,
                    help=( "comma-separated features to ignore (i.e. will be ignored; will not be collased as 'others')."
                           " default='chromosome,protein'" ),
                    default="chromosome,protein")
parser.add_argument("--ignore-intergenic", action="store_false", default=False, target="report_intergenic",
                    help="don't report intergenic features (i.e. features that intersect nothing or only 'chromosome')")
args = parser.parse_args()

## get new priorities
new_priority = ([] if args.priority is None else args.priority.split(','))
new_priority_last = ([] if args.priority_last is None else [feature for feature in args.priority_last.split(',') if feature not in new_priority])
priority_order = [feature for feature in priority_order if (feature not in new_priority and feature not in new_priority_last)]
priority_last = [feature for feature in priority_last if (feature not in new_priority and feature not in new_priority_last)]
priority_order = new_priority + priority_order
priority_last = priority_last + new_priority_last

## get new ignore
ignore = (ignore if args.ignore is None else args.ignore.split(','))
ignore_intergenic = args.ignore_intergenic

## see if we can use pybedtools to standardise output so we don't have to run intersect for gff and bed separately
