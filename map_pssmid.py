import re
import sys
sys.path.append("/mnt/chaelab/rachelle/src")
from data_manip import *

cdd_versions_fname = "/mnt/chaelab/rachelle/data/cdd/cdd.versions.tsv"
## cdd_versions_fname should point to a tsv with headers (important columns are: pssmid, shortn)

def map_pssmid(fname, pssm_col, fout, header = False, output_header = False, sep = '\t',
               new_col = "name", cdd_versions_fname = cdd_versions_fname, cdd_tbl_fname = None,
               pattern = "\d+", remove_oldcol = True):
    data = [x.split(sep) for x in splitlines(fname) if x]
    # with open(fname, 'r') as f:
    #     data = [x[:-1].split(sep) for x in f.readlines() if x]
    ## parse header (if present)
    if header and \
       re.search(pattern, data[0][pssm_col]) and \
       re.search(pattern, data[0][pssm_col]).group(0).isnumeric():
        print("Number detected in header row, PSSM-Id column. Parsing first row as data.")
        header = False
    elif header:
        data_header = data[0]
        data = data[1:]
        if remove_oldcol:
            header_to_keep = [i for i, v in enumerate(data_header) if v != new_col]
            data = [[line[i] for i in header_to_keep] for line in data]
            data_header = [data_header[i] for i in header_to_keep]
    ## read cdd_versions/tbl
    if cdd_tbl_fname is not None:
        with open(cdd_tbl_fname, 'r') as f:
            cdd_data = [x[:-1].split('\t') for x in f.readlines()]
            ## columns: pssmid, domain id, shortname, long name, length
            cdd_data = {int(x[0]): x[2] for x in cdd_data}
    else:
        with open(cdd_versions_fname, 'r') as f:
            cdd_data = [x[:-1].split('\t') for x in f.readlines() if x]
            cdd_get = make_custom_get(cdd_data[0])
            cdd_data = {cdd_get(x, "pssmid"): cdd_get(x, "shortname") for x in cdd_data[1:]}
    ## map shortname to pssmid
    data = [x + [cdd_data.get(int(re.search(pattern, x[pssm_col]).group(0)), "NA")] for x in data]
    ## write out
    if header and output_header:
        data = [data_header + [new_col]] + data
    with open(fout, "w+") as f:
        f.write('\n'.join(['\t'.join([str(y) for y in x]) for x in data]))
    return
