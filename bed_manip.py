import csv
import subprocess

def bed_to_list(bed_f, row_start = 0):
    with open(bed_f, 'r') as f:
        raw = [x[:-1].split('\t') for x in f.readlines()[row_start:]]
        bed_data = list(map(lambda x: [x[0], x[1] if not x[1].isdigit() else int(x[1]),
                                       x[2] if not x[2].isdigit() else int(x[2])] + x[3:],
                            filter(lambda x: len(x) > 1, raw)))
    return bed_data

# extends bed region by ext_len in both directions
# lower bound is 0 and upper bound is chr size
def extend(fname, fout, fchr, ext_len=20000):
    """
    :fname: path to bed file with regions to be extended
    :fout: output file path
    :fchr: path to bed file with chromosomes only
    :ext_len: extension size (in bases)
    """
    chr_sizes = {}
    with open(fchr, "r") as f:
        lines = [x[:-1].split("\t") for x in f.readlines()]
        for line in lines:
            chr_sizes[line[0]] = int(line[2])
        
    output = [] # list to store all output lines
    with open(fname, "r") as f:
        lines = [x[:-1].split("\t") for x in f.readlines()]
        print(len(lines))
        for line in lines:
            line[1] = str(max(0,
                              int(line[1])-ext_len))
            line[2] = str(min(chr_sizes[line[0]],
                              int(line[2])+ext_len))
            output.append('\t'.join(line))

    # write features to file
    import os
    f = open(fout, "w+")
    f.write('\n'.join(output))
    f.close()
    # sort for thoroughness
    os.system("bedtools sort -i {} > temp.txt; mv temp.txt {}".format(fout, fout))
    return 0

def extract_bed_generic(genes, bed, out_dir, feature = "gene", remove_chr = True, add_chr = "", encoding = "utf-8", for_vcftools = False):
    print("Extracting gene entries of type " + feature + " from " + bed)
    ## grep and awk relevant genes & feature entries
    grep = subprocess.Popen(("grep", "-E", "|".join(genes), bed), stdout=subprocess.PIPE)
    bed_raw = [entry.split('\t') for entry in
               subprocess.check_output(("awk", "$8 == \"" + feature + "\""),
                                       stdin=grep.stdout).decode(encoding).split('\n') if entry]
    ## modify chromosome column as indicated by arguments
    bed_raw = [[add_chr + x[0] if not remove_chr else re.search("[Cc]hr(.+)",x[0]).group(1)] + x[1:]
               for x in bed_raw]
    ## write reduced bed file
    tmp_bed = out_dir + "/ranges.bed"
    ## header is required by vcftools although it is not necessary for regular BED files
    list_to_bed(bed_raw, tmp_bed, for_vcftools = for_vcftools)
    ## return location of reduced bed file
    print("Successfully extracted gene feature positions into " + tmp_bed)
    if for_vcftools:
        print("(Note that this bed file is unconventional in that it has a header that is required by vcftools)")
    return tmp_bed

def extract_bed(genes, bed, out_dir, **for_extract_bed_generic):
    return extract_bed_generic(genes, bed, out_dir, for_vcftools = False, **for_extract_bed_generic)

def extract_bed_for_vcftools(genes, bed, out_dir, **for_extract_bed_generic):
    return extract_bed_generic(genes, bed, out_dir, for_vcftools = True, **for_extract_bed_generic)

def list_to_tsv(l, fname):
    f = open(fname, "w+")
    f.write('\n'.join(['\t'.join(str(x)) for x in l]))
    f.close()

def list_to_bed(l, fname, for_vcftools = False):
    if for_vcftools:
        l = [["chrom", "chromStart", "chromEnd"]] + l
    list_to_tsv(l, fname)
