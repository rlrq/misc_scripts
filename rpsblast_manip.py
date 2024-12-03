# import sys
# sys.path.append("/mnt/chaelab/rachelle/src")

import itertools

def index_cdd_versions(fname):
    """
    Create dictionary of {'<pssmid>': '<shortname>'} for live PSSMs from cdd.versions file
    """
    ## get length of columns
    header = ''
    with open(fname, 'r') as f:
        for entry in f:
            if entry[0] == "#":
                header = entry[:-1]
                break
    shortname_start = header.index("ShortName")
    pssmid_start = header.index("PssmId")
    lv_start = header.index("Lv")
    shortname_end = pssmid_start - 1
    pssmid_end = header.index("Parent") - 1
    lv_end = header.index("Rl") - 1
    def get_shortname(entry):
        return entry[shortname_start:shortname_end].rstrip()
    def get_pssmid(entry):
        return entry[pssmid_start:pssmid_end].rstrip()
    def get_lv(entry):
        return entry[lv_start:lv_end].rstrip()
    ## index shortnames by pssmid
    output = {}
    with open(fname, 'r') as f:
        for entry in f:
            if get_lv(entry) == '1':
                output[get_pssmid(entry)] = get_shortname(entry)
    return output

def index_cdd_tbl(fname):
    """
    Create dictionary of {'<pssmid>': '<shortname>'} for live PSSMs from cddid.tbl file
    """
    output = {}
    with open(fname, 'r') as f:
        for entry in f:
            output[entry.split('\t')[0]] = entry.split('\t')[2]
    return output

def get_pssmid(sseqid):
    return sseqid.split(':')[-1].split('|')[-1]

def filter_rpsblast(blasttsv, fout, pssm_ids, col_names = None,
                    cols_to_write = None, col_names_out = None,
                    write_domain_name = False, cdd_versions = None, cdd_tbl = None) -> None:
    """
    Filter RPS-BLAST output for specific PSSM-IDs and write to file.
    
    Arguments:
        blasttsv (str): path to RPS-BLAST output tsv file (outfmt 6)
        fout (str): path to output file
        pssm_ids (list/tuple/set): PSSM-Ids to keep. If empty, keeps all.
        col_names (list/tuple): column names if blasttsv has non-standard output field combo/order.
            Required columns: sseqid.
            If not provided, assumes default:
            'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
        cols_to_write (list/tuple): optional, indices of columns to write to file.
            If not provided, all columns will be kept.
        col_names_out (list/tuple): optional, header (colum names) in output file.
            If not provided, no header will be written.
        write_domain_name (bool): default=False. If True, writes PSSM-Id's short name in last column.
            Requires cdd_versions OR cdd_tbl if True.
        cdd_versions (str): path to cdd.versions file. Used to retrieve short name for PSSM-Ids.
            Required if write_domain_name=True and cdd_tbl=None.
        cdd_tbl (str): path to cddid.tbl file. Used to retrieve short name for PSSM-Ids.
            Required if write_domain_name=True and cdd_versions=None.
    """
    pssm_ids = set(map(str, pssm_ids))
    ## index pssm ids
    if write_domain_name:
        if cdd_tbl:
            indexed_shortnames = index_cdd_tbl(cdd_tbl)
        elif cdd_versions:
            indexed_shortnames = index_cdd_versions(cdd_versions)
        else:
            print("cdd_tbl or cdd_versions required if write_domain_name=True")
            return
        get_shortname = lambda pssmid: [indexed_shortnames[pssmid]]
    else:
        get_shortname = lambda pssmid: []
    ## index-related stuff
    pssmid_i = 1 if not col_names else col_names.index("sseqid")
    if cols_to_write:
        filter_cols = lambda entry, pssmid: [entry[i] for i in cols_to_write] + get_shortname(pssmid)
    else:
        filter_cols = lambda entry: entry + get_shortname(pssmid)
    ## parse data
    output = []
    with open(blasttsv, 'r') as f:
        for entry in f:
            entry = entry[:-1].split('\t')
            pssmid = get_pssmid(entry[pssmid_i])
            if pssm_ids and (pssmid not in pssm_ids):
                continue
            else:
                output.append(filter_cols(entry, pssmid))
    ## write
    if col_names_out:
        output = [col_names_out] + output
    with open(fout, 'w+') as f:
        for line in output:
            f.write('\t'.join(line) + '\n')
    return

def rps_to_nlr(blasttsv, fout, code, col_names = None,
               by_start = True, collapse = True, collapse_order = None,
               merge_overlapping_with_same_code = True):
    """
    Filter RPS-BLAST output for specific PSSM-IDs and write to file.
    
    Arguments:
        blasttsv (str): path to RPS-BLAST output tsv file (outfmt 6)
        fout (str): path to output file
        code (dict): dict of
            {'<short representation of domain; single letter recommended>': [<pssmid1>, <pssmid2>]}.
            (e.g. {'T': [366714], 'C': [271353, 375519], 'R': [368547], 'N': [366375], 'L': [227223]})
        col_names (list/tuple): column names if blasttsv has non-standard output field combo/order.
            Required columns: qseqid, sseqid, qstart, qend.
            If not provided, assumes default:
            'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
        by_start (bool): write column of domain architecture in order of start position
        collapse (bool): write column of collapsed domain architecture. Repeated domains will only be
            represented by their code once. (e.g. TNNL --> TNL)
        collapse_order (tuple/list): order in which to write collapsed architecture.
            Used only if collapse=True. If not provided, defaults to alphabetical sorting.
            (e.g. if collapse_order = ['T', 'R', 'C', 'N', 'L'], TNTNLN --> TNL)
        merge_overlapping_with_same_code (bool): if 2 PSSMs of same code overlap, treat as single domain.
            Used only if by_start = True.
            (e.g. given {'A': [pssmid1, pssmid2], 'B': [pssmid3]} 
            and that pssmid1 spans 5aa-10aa and pssmid2 spans 6aa-10aa and pssid3 spans 15aa-20aa,
            domain architecture is AB if merge_overlaping_with_same_code=True and 
            AAB if merge_overlapping_with_same_code=False.)
    """
    ## parse pssmids to str
    code = {k: tuple(map(str, v)) for k, v in code.items()}
    pssmids = set(itertools.chain(*code.values()))
    inv_code = {}
    for k, v in code.items():
        for pssmid in v:
            inv_code[pssmid] = k
    ## index-related stuff
    pssmid_i = 1 if not col_names else col_names.index("sseqid")
    qseqid_i = 0 if not col_names else col_names.index("qseqid")
    qstart_i = 6 if not col_names else col_names.index("qstart")
    qend_i = 7 if not col_names else col_names.index("qend")
    ## filter for desired pssmids and group by qseqid
    dat = {}
    with open(blasttsv, 'r') as f:
        for line in f:
            entry = line[:-1].split('\t')
            pssmid = get_pssmid(entry[pssmid_i])
            if pssmid in pssmids:
                dat[entry[qseqid_i]] = dat.get(entry[qseqid_i], []) + \
                                       [[int(entry[qstart_i]), int(entry[qend_i]),
                                         inv_code[pssmid]]]
    ## sort and collapse
    output = {}
    start_i, end_i, code_i = range(3)
    cols_to_write = []
    ## by_start
    if by_start:
        col_name = "by_start"
        cols_to_write.append(col_name)
        d_by_start = {}
        for qseqid, qseq_dat in dat.items():
            qseq_dat.sort() ## default sort order as entry elements are already in sorting order
            architecture = ''
            last_end = -1
            last_code = ''
            for entry in qseq_dat:
                if (int(entry[start_i]) > int(last_end)
                    or not merge_overlapping_with_same_code
                    or (merge_overlapping_with_same_code and entry[code_i] != last_code)):
                    architecture += entry[code_i]
                last_end = entry[end_i]
                last_code = entry[code_i]
            d_by_start[qseqid] = architecture
        output[col_name] = d_by_start
    ## collapse
    if collapse:
        col_name = "collapse"
        cols_to_write.append(col_name)
        d_collapse = {}
        sort_key = (lambda x: x) if not collapse_order else (lambda x: collapse_order.index(x))
        for qseqid, qseq_dat in dat.items():
            qseq_codes = list(set(entry[code_i] for entry in qseq_dat))
            qseq_codes.sort(key = sort_key)
            d_collapse[qseqid] = ''.join(qseq_codes)
        output[col_name] = d_collapse
    ## write
    with open(fout, "w+") as f:
        qseqids = tuple(sorted(dat.keys()))
        f.write('\t'.join(["qseqid"] + cols_to_write) + '\n')
        for qseqid in qseqids:
            f.write('\t'.join([qseqid] + [output[col_name][qseqid] for col_name in cols_to_write]) + '\n')
    # return output
    return
