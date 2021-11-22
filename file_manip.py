import re
import sys
sys.path.append("/mnt/chaelab/rachelle/src")
import fileinput
from data_manip import splitlines
from fasta_manip import *

def filter_column(fname, fitems, col, fout, inverse = False, delim = '\t'):
    with open(fname, 'r') as f:
        dat = [x[:-1].split(delim) for x in f.readlines()]
    items = []
    with open(fitems, 'r') as f:
        items_raw = [x[:-1].split(delim) for x in f.readlines()]
        for item in items_raw:
            items.extend(item)
    output = [delim.join(x) for x in dat if x[col] in items] if not inverse else \
             [delim.join(x) for x in dat if x[col] not in items]
    with open(fout, "w+") as f:
        f.write('\n'.join(output))
    return

## set exists to False to create file if it doesn't exist; else, error raised
def prepend(fname, to_prepend, exists = True):
    ## read
    try:
        with open(fname, 'r') as original:
            data = original.read()
    except FileNotFoundError as e:
        if exists: raise e
        else: data = ''
    ## write
    with open(fname, 'w+') as modified:
        modified.write(to_prepend + data)
    return
    
# fname = "/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_gene.bed"
# fitems = "/mnt/chaelab/shared/NLR/geneID/all_NLR.txt"
# col = 3
# fout = "/mnt/chaelab/shared/NLR/bed/all_NLR.bed"
# ## filter_column(fname, fitems, col, fout)

# fasta = "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta"
# fout = "/mnt/chaelab/shared/NLR/sequence/all_NLR.fasta"
def intersect_bed_fasta(bed, fasta, fout):
    with open(bed, 'r') as f:
        bed_dat = [x[:-1].split('\t') for x in f.readlines()]
    fasta_dat = fasta_to_dict(fasta)
    output = {}
    for feat in bed_dat:
        chrom = feat[0]
        start, end = feat[1:3]
        geneID = feat[3]
        feat_name = geneID + '_' + chrom + ':' + start + '-' + end
        output[feat_name] = fasta_dat[chrom][int(start):int(end)]
    dict_to_fasta(output, fout)


def single_col_to_list(fname, delim = '\t'):
    """
    reads a single-column file and generates a list of each line (removes trailing newline char)
    """
    raw = splitlines(fname)
    return [x for x in raw]

def cat_files(fnames, fout):
    with open(fout, "w+") as f, fileinput.input(fnames) as fin:
        for line in fin:
            f.write(line)
    return fout
    
# fname = "/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70_AT5G58120.Col-0.blastn_summary.tsv"
# fout = "/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70_AT5G58120.Col-0.blastn_topHits.tsv"
# bs_col = -1
# keep_best_hit_per_query_subject_pair(fname, fout, -1, min_length = 500)
def keep_best_hit_per_query_subject_pair(fname, fout, bs_col, q_col = 0, s_col = 1, delim = '\t', header = True, length_col = 3, min_length = 0):
    with open(fname, 'r') as f:
        data_raw = [x[:-1].split(delim) for x in f.readlines()]
    output = [delim.join(data_raw[0])] if header else []
    data = data_raw[1:] if header else data_raw
    data_rows = len(data)
    data.sort()
    last_row = data[0]
    last_q, last_s, last_bs = last_row[q_col], last_row[s_col], float(last_row[bs_col])
    for i,row in enumerate(data[1:]):
        if int(row[length_col]) < min_length:
            continue
        curr_q, curr_s, curr_bs = row[q_col], row[s_col], float(row[bs_col])
        if curr_q != last_q or curr_s != last_s:
            output.append(delim.join(last_row))
            last_row = row
            last_q, last_s, last_bs = last_row[q_col], last_row[s_col], float(last_row[bs_col])
        else:
            if curr_bs > last_bs:
                last_row = row
                last_q, last_s, last_bs = last_row[q_col], last_row[s_col], float(last_row[bs_col])
        if i == data_rows - 2:
            output.append(delim.join(last_row))
    f = open(fout, "w+")
    f.write('\n'.join(output))
    f.close()

# fname = "/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70_AT5G58120.Col-0.blastn_topHits.tsv"
# fout = "/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70_AT5G58120.Col-0.blastn_topHits_AT5G58120.tsv"
# subjects = ["AT5G58120_Chr5:23517491-23521067"]
# keep_query_aligned_to_selected_subjects(fname, fout, subjects)
def keep_query_aligned_to_selected_subjects(fname, fout, subjects, q_col = 0, s_col = 1, delim = '\t', header = True):
    with open(fname, 'r') as f:
        data_raw = [x[:-1].split(delim) for x in f.readlines()]
    output = [delim.join(data_raw[0])] if header else []
    data = data_raw[1:] if header else data_raw
    selected_queries = []
    for subject in subjects:
        selected_queries.extend([x[q_col] for x in data if subject in x[s_col]])
    output.extend([delim.join(x) for x in data if x[q_col] in selected_queries])
    f = open(fout, "w+")
    f.write('\n'.join(output))
    f.close()

def replace_ext(fname, new_ext, old_ext = None):
    if not old_ext:
        fname = re.search("^.+(?=\.[^.]+$)", fname).group(0)
    else:
        fname = re.search(f"^.+(?=\.{old_ext}$)", fname).group(0)
    return f"{fname}.{new_ext}"
