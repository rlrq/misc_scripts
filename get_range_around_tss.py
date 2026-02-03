#!/usr/bin/python3
import argparse
import re

import sys
sys.path.append("/mnt/chaelab/rachelle/src")
import range_manip as rm

parser = argparse.ArgumentParser(description="get range around TSS(s) of genes")
parser.add_argument("gff", type=str, help="path to GFF3 file")
parser.add_argument("fout", type=str, help="path to output BED file")
parser.add_argument("--separate", help="do not merge TSS ranges of different mRNA of same gene",
                    action="store_true", default=False)
parser.add_argument("--up", type=int, help="upstream range extension (bp)", default=2000)
parser.add_argument("--down", type=int, help="downstream range extension (bp)", default=2000)
args = parser.parse_args()

fout = args.fout
up_range = args.up
down_range = args.down

print("Reading GFF3 file")
with open("/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.gff", 'r') as f:
# with open(args.gff, 'r') as f:
    ## get mRNA entries and group by gene
    mrna_entries = {}
    for entry in f:
        ## skip if header
        if entry[0] == '#': continue
        ## skip if not mRNA
        split_line = entry.rstrip().split('\t')
        if split_line[2] != "mRNA": continue
        ## store by gene ID
        gene = re.search("(?<=^Parent=)[^;]+|(?<=;Parent=)[^;]+", split_line[8])[0]
        mrna_entries[gene] = mrna_entries.get(gene,[]) + [tuple(split_line)]

## get TSS of each mRNA
print("Getting TSS range of each mRNA")
tss_ranges = {}
for gene, entries in mrna_entries.items():
    tss_ranges_gene = tss_ranges.get(gene, {})
    for entry in entries:
        mrna_id = re.search("(?<=^ID=)[^;]+|(?<=;ID=)[^;]+", entry[8])[0]
        chrom = entry[0]
        strand = entry[6]
        if strand == '+':
            tss_pos = int(entry[3])
            tss_ranges_gene[mrna_id] = (tss_pos-up_range, tss_pos+down_range)
        elif strand == '-':
            tss_pos = int(entry[4])
            tss_ranges_gene[mrna_id] = (tss_pos-down_range, tss_pos+up_range)
    tss_ranges[gene] = tss_ranges_gene

##### WRITE
print("Writing ranges to file")
## merge TSS extended ranges if not --separate
i = 0
if not args.separate:
    with open(fout, 'w+') as f:
        for gene, mrna_dat in tss_ranges.items():
            if i % 1000 == 0: print(i)
            i += 1
            base_entry = mrna_entries[gene][0]
            chrom = base_entry[0]
            source = base_entry[1]
            strand = base_entry[6]
            merged_ranges = rm.ranges_union([list(mrna_dat.values())])
            ## range should be converted to 0-index for BED output
            for r in merged_ranges:
                _ = f.write('\t'.join(map(str, [chrom, max(0, r[0]-1), r[1], f"{gene}_TSS_region",
                                                '.', strand, source, "TSS_region", '.',
                                                f"Gene={gene}"])) + '\n')
else:
    with open(fout, 'w+') as f:
        for gene, mrna_dat in tss_ranges.items():
            if i % 1000 == 0: print(i)
            i += 1
            base_entry = mrna_entries[gene][0]
            chrom = base_entry[0]
            source = base_entry[1]
            strand = base_entry[6]
            ## range should be converted to 0-index for BED output
            for mrna_id, tss_range in mrna_dat.items():
                _ = f.write('\t'.join(map(str, [chrom, max(0, tss_range[0]-1), tss_range[1],
                                                f"{mrna_id}_TSS_region",
                                                '.', strand, source, "TSS_region", '.',
                                                f"mRNA={mrna_id};Gene={gene}"])) + '\n')
