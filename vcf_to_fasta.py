#!/usr/bin/python3
import argparse

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import vcf

from Bio import Seq

from annotation import GFF
from index_fasta import IndexedFasta
from fasta_manip import dict_to_fasta

## function to check if position is in (sorted) chromosome ranges
def pos_in_range(chrom, pos, ranges):
    if chrom not in ranges: return False
    for r in ranges[chrom]:
        if pos < r[0]: return False
        elif pos >= r[0] and pos < r[1]: return True
    return False

## return indexed output ranges + strands
def parse_input_ranges(f_bed, adjust_dir = False):
    ## get ranges to output
    output_ranges = GFF(fname = f_bed, fmt = "BED")
    output_ranges_indexed = {}
    strands = {}
    for entry in output_ranges:
        strands[entry.molecule] = strands.get(entry.molecule, set()).union({entry.strand})
        output_ranges_indexed[entry.molecule] = (output_ranges_indexed.get(entry.molecule, []) +
                                                 [(entry.start, entry.end)])
    ## if adjust_dir = True, check if all entries are of the same strand. fail otherwise.
    if adjust_dir:
        for chrom, strand_dat in strands.items():
            if len(strand_dat) >= 2:
                raise Exception(("All BED entries must be of the same strand if adjust_dir=True."
                                 f" Check molecule {chrom}."))
    ## sort indexed output ranges
    output_ranges_indexed = {k: sorted(v) for k, v in output_ranges_indexed.items()}
    return output_ranges_indexed, strands

## main function
def vcf_to_seq(f_bed, f_fasta, f_vcf, f_out, feature_id = None,
               adjust_dir = True, remove_redundant = False, homozygous = False):
    output_ranges_indexed, strands = parse_input_ranges(f_bed, adjust_dir = adjust_dir)
    ## function to check if pos is in output ranges
    pos_in_output = lambda chrom, pos: pos_in_range(chrom, pos, output_ranges_indexed)
    ## parse & filter variants
    variants = list(record for record in vcf.Reader(filename = f_vcf) if pos_in_output(record.CHROM, record.POS))
    ## parse FASTA file
    ifasta = IndexedFasta(f_fasta)
    ## get reference sequence for each bed entry
    seq_ref = {}
    for chrom, chrom_ranges in output_ranges_indexed.items():
        seq_ref[chrom] = {}
        for chrom_range in chrom_ranges:
            seq_ref[chrom][chrom_range] = ifasta[chrom][chrom_range[0]:chrom_range[1]]
    ## create variant sequences for each entry
    seq_variants = {}
    for chrom, chrom_data in seq_ref.items():
        seq_variants[chrom] = {}
        for chrom_range, seq in chrom_data.items():
            start, end = chrom_range
            seq_variants[chrom][chrom_range] = {}
            curr_data = seq_variants[chrom][chrom_range]
            ## get reference sequence
            ref_seq = seq_ref[chrom][chrom_range]
            ## filter for variants in range and sort by position
            variants_in_range = sorted((var for var in variants if
                                        (var.CHROM == chrom and var.POS >= start and var.POS < end)),
                                       key = lambda var: var.POS)
            ## iterate through all variants in range
            for variant in variants_in_range:
                offset_pos = variant.POS - start
                for sample in variant.samples:
                    sample_name = sample.sample
                    ## add reference sequence to output if sample not encountered yet
                    if sample_name not in curr_data:
                        if homozygous:
                            curr_data[sample_name] = [list(ref_seq)]
                        else:
                            ## store two sequences to accommodate heterozygosity, assume phased
                            curr_data[sample_name] = [list(ref_seq), list(ref_seq)]
                    ## replace reference sequence w/ called base
                    try:
                        called_bases = sample.gt_bases.split('/')
                    except:
                        called_bases = ['N', 'N']
                    ## first called base is processed the same whether homozygosity is assumed
                    curr_data[sample_name][0][offset_pos] = called_bases[0]
                    ## second called base only processed when homozygosity is not assumed
                    if not homozygous:
                        curr_data[sample_name][1][offset_pos] = called_bases[1]
    ## get all sample names
    samples = set().union(*(range_data.keys()
                            for chrom_data in seq_variants.values()
                            for range_data in chrom_data.values()))
    ## reconstitute sequences for all samples
    seq_pseudomolecule = {}
    for chrom, chrom_ranges in output_ranges_indexed.items():
        seqs_curr = {}
        for chrom_range in sorted(chrom_ranges):
            ref_seq = [seq_ref[chrom][chrom_range],seq_ref[chrom][chrom_range]] ## accommodate het
            seq_samples = seq_variants.get(chrom, {}).get(chrom_range, {})
            for sample in samples:
                if homozygous:
                    seqs_curr[f"{sample}|0"] = (seqs_curr.get(f"{sample}|0", '') +
                                                ''.join(seq_samples.get(sample, ref_seq)[0]))
                else:
                    ## accommodate heterozygotes
                    for i in (0,1):
                        seqs_curr[f"{sample}|{i}"] = (seqs_curr.get(f"{sample}|{i}", '') +
                                                      ''.join(seq_samples.get(sample, ref_seq)[i]))
        ## remove redundant
        if remove_redundant and not homozygous:
            for sample in samples:
                if seqs_curr[f"{sample}|0"] == seqs_curr[f"{sample}|1"]:
                    del seqs_curr[f"{sample}|1"]
        ## add id/ranges to seq name
        if feature_id is not None:
            seqs_curr = {f"{k}|{feature_id}": seq for k, seq in seqs_curr.items()}
        else:
            str_ranges = ','.join('-'.join(map(str, r)) for r in chrom_ranges)
            seqs_curr = {f"{k}|{chrom}:{str_ranges}": seq for k, seq in seqs_curr.items()}
        ## adjust direction
        if adjust_dir:
            strand = strands[chrom]
            if strand == '-':
                seqs_curr = {f"{k}|revcomp": Seq.Seq(seq).reverse_complement() for k, seq in seqs_curr.items()}
        ## add to output
        seq_pseudomolecule = {**seq_pseudomolecule, **seqs_curr}
    ## write to fasta file
    dict_to_fasta(seq_pseudomolecule, f_out)
    return


parser = argparse.ArgumentParser(description="create pseudomolecules from VCF and FASTA files")
parser.add_argument("--vcf", type=str, help="path to VCF file")
parser.add_argument("--fasta", type=str, help="path to FASTA file")
parser.add_argument("--bed", type=str, help="path to BED file containing regions to subset")
parser.add_argument("--fout", type=str, help="output file")
parser.add_argument("--feature-id", type=str, help="ID of regions, used to name output sequences; subset coordinates used if not provided", default=None, dest="feature_id")
parser.add_argument("--adj-dir", "--adjust-dir", help="output reverse complement if on minus strand", action="store_true", default=False, dest="adjust_dir")
parser.add_argument("--rm-redundant", "--remove-redundant", help="remove redundant sequences (e.g. when individual is homozygous at all sites, only one sequence will be retained)", action="store_true", default=False, dest="remove_redundant")
parser.add_argument("--homo", "--homozygous", help="assume all individuals are homozygous (i.e. only parses first reported base call at each position, outputs one sequence per individual)", action="store_true", default=False, dest="homozygous")
args = parser.parse_args()

vcf_to_seq(args.bed, args.fasta, args.vcf, args.fout, feature_id = args.feature_id,
           adjust_dir = args.adjust_dir, remove_redundant = args.remove_redundant,
           homozygous = args.homozygous)


# f_fasta = "/mnt/chaelab/rachelle/data/Samericanum/lin2023/sp1102.v5.fa"
# f_bed = "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/bed/PAP2_sp1102.final.CDS.bed"
# f_vcf = "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/vcf/PAP2_sp1102.final.CDS.lin2023.recode.vcf"
# f_out = "/mnt/chaelab/rachelle/zzOtherzz/DM3/2023_manuscript/seq/PAP2.lin2023.CDS.fasta"
# feature_id = "sp1102chr01.1927.t1_sp1102chr01.1931.t1" ## used only to name output sequences

# adjust_dir = True
# remove_redundant = True
