import os
import re
import subprocess
from basics import *
from data_manip import *

def mummer_intersect(mums, bed, fout):
    
    ## parse bed file
    ## (retain only first 3 fields for sole purpose of eliminating non-overlapping mums entries)
    ## (actual bedtools intersect will be carried out on original file)
    print("Parsing BED file")
    bed_dict = {}
    for entry in [line.split('\t') for line in splitlines(bed)]:
        chrom = entry[0]
        bed_dict[chrom] = bed_dict.get(chrom, []) + [list(map(int, entry[1:3]))]
    bed_dict = {k: sorted(v, key = lambda x: int(x[1])) for k, v in bed_dict.items()}
    
    ## parse mums file
    print("Parsing MUMS file")
    mums_bed = []
    nonref = ''
    forward = True
    i = 0
    for line in open(mums, 'r').readlines():
        i += 1
        if i % 100000 == 0:
            print(f"Processing line {i}")
        entry = [field for field in ''.join([c for c in line[2:] if c != '\n']).split(' ') if field]
        ## if is header (indicating molecule in non-reference genome)
        if line[0] == '>':
            nonref = entry[0]
            print(f"Non-reference molecule: {nonref}")
            forward = len(entry) == 1
            continue
        else:
            refchrom = entry[0]
            refstart, nonrefstart, length = map(int, entry[1:])
            if any_overlap((refstart - 1, refstart + length - 1), bed_dict.get(refchrom, tuple())):
                mums_bed.append([refchrom, refstart - 1, refstart + length - 1,
                                 nonref, nonrefstart - 1, nonrefstart + length - 1,
                                 '+' if forward else '-'])
                
    ## bedtools intersect
    print("Executing bedtools intersect")
    tmp_mums_bed = fout + '.tmp'
    open(tmp_mums_bed, "w+").write('\n'.join(['\t'.join(map(str, line)) for line in mums_bed]))
    fout_f = open(fout, "w+")
    subprocess.run(("bedtools", "intersect", "-wao", "-a", bed, "-b", tmp_mums_bed), stdout = fout_f)
    fout_f.close()
    os.remove(tmp_mums_bed)
    return
