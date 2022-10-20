import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from fasta_manip import fasta_to_dict, dict_to_fasta
from data_manip import splitlines

def blast6_to_qseqid_dict(fname, *cols_to_keep, i_qseqid = 0) -> dict:
    """
    Read BLAST output file (fmt 6; tsv) into dictionary indexed by
    qseqid. Format: {<qseqid>: [[<hit1 data>], [<hit2 data>], ...]}
    
    Arguments:
        fname (str): path to BLAST fmt 6 (tsv) output file
        *cols_to_keep (int): indices of columns to keep; data will
            be ordered in the same order as specified here; a field
            may be specified multiple times
        i_qseqid (int): column index of qseqid field
            (default=0; per fmt 6 default field)
    
    Returns:
        dict
    """
    ## make keep_cols function to generate retained data
    if cols_to_keep:
        def keep_cols(line_dat):
            return [line_dat[i] for i in cols_to_keep]
    else:
        def keep_cols(line_dat):
            return line_dat
    ## parse data
    dat = {}
    with open(fname, 'r') as f:
        for line in f:
            line_dat = line[:-1].split('\t')
            qseqid = line_dat[i_qseqid]
            dat[qseqid] = dat.get(qseqid, []) + [keep_cols(line_dat)]
    return dat

def blast6_top_n(blast, n = 5, fout = None, fields = None) -> None:
    """
    Retain only top 'n' hit for each qseqid (sorted by bitscore). If
    tie for Nth place, all hits with threshold score are retained.
    
    Arguments:
        blast (str): path to BLAST fmt 6 (tsv) output file
        n (int): top N number of hits to keep per qseqid (default=5)
        fout (str): optional, path to output file
        fields (list/tuple): rblast field names, required only if
            non-standard order/combo of fields in rblast
    """
    if fields:
        i_qseqid = fields.index("qseqid")
        i_bitscore = fields.index("bitscore")
    else:
        i_qseqid = 0
        i_bitscore = -1
    dat = blast6_to_qseqid_dict(blast, i_qseqid = i_qseqid)
    for qseqid, qdat in dat.items():
        if len(qdat < n):
            continue
        threshold_score = float(sorted(qdat, key = lambda x: float(x[i_bitscore]))[-n][i_bitscore])
        dat[qseqid] = [x for x in qdat if float(x[i_bitscore]) >= threshold_score]
    return dat

def blast6_top_n_memsave(blast, fout = None, n = 5, fields = None) -> dict:
    """
    Retain only top 'n' hit for each qseqid (sorted by bitscore). If
    tie for Nth place, all hits with threshold score are retained.
    
    Memory saving version; may take longer to run as threshld score is
    dynamically calculated after every hit is processed but should save
    memory by not reading all reads to memory first before filtering.
    
    Arguments:
        blast (str): path to BLAST fmt 6 (tsv) output file
        fout (str): optional, path to output file
        n (int): top N number of hits to keep per qseqid
        fields (list/tuple): rblast field names, required only if
            non-standard order/combo of fields in rblast
    
    Returns:
        dict
    """
    ## parse fields
    if fields:
        i_qseqid = fields.index("qseqid")
        i_bitscore = fields.index("bitscore")
    else:
        i_qseqid = 0
        i_bitscore = -1
    ## parse data
    dat = {}
    threshold = {}
    with open(blast, 'r') as f:
        for line in f:
            line_dat = line[:-1].split('\t')
            qseqid = line_dat[i_qseqid]
            ## if no entry for this qseqid, make it
            if qseqid not in dat:
                dat[qseqid] = [line_dat]
                threshold[qseqid] = -1 ## start with nonexistent threshold for seeding
            ## if entry for this qseqid already exists, compare w/ existing threshold & data
            else:
                bitscore = float(line_dat[i_bitscore])
                ## skip if bitscore lower than threshold
                if bitscore < threshold[qseqid]:
                    continue
                ## add if bitscore equal to threshold or not yet 'n' number of hits
                elif (bitscore == threshold[qseqid] or
                      len(dat[qseqid]) < n):
                    dat[qseqid].append(line_dat)
                ## add and update threshold if bitscore > threshold and already 'n' number of hits
                else:
                    ## update threshold
                    new_threshold = sorted([float(x[i_bitscore]) for x in (dat[qseqid] + [line_dat])])[-n]
                    threshold[qseqid] = new_threshold
                    ## filter by threshold
                    dat[qseqid] = [x for x in (dat[qseqid] + [line_dat]) if float(x[i_bitscore]) >= new_threshold]
    ## write
    if fout:
        with open(fout, "w+") as f:
            for qseqid, qdat in dat.items():
                for hit_dat in qdat:
                    f.write('\t'.join(hit_dat) + '\n')
    return dat

## required fields: qseqid, sseqid, bitscore
def filter_from_blast6(rblast, fout, fields = None, targets = None,
                       eq = lambda qseqid, sseqid: qseqid == sseqid) -> None:
    """
    Retain only top hit for each qseiqd (sorted by bitscore).
    Outputs a tsv file mapping qseqid to best matching sseqid + bitscore.
    If 'targets' is specified, retain only qseqid which best matching
    hit is found in 'targets'.
    
    Arguments:
        rblast (str): path to BLAST fmt 6 (tsv) output file
        fout (str): path to output file.
            Fmt: <qseqid>\t<sseqid of best hit>\t<bitscore>[\t<sseq is target (if targets!=None)>]
        fields (list/tuple): rblast field names, required only if
            non-standard order/combo of fields in rblast
        targets (list/tuple/set): optional, IDs of target genes/features
        eq (func): optional, function to equate qseqid and sseqid
            if they follow different formats.
            E.g. If qseqids are isoform names (e.g. 'AT1G10010.1')
            and sseqids are gene IDs with extraneous fields
            (e.g. 'Ref|AT1G10010|CDS') you'd need to use this variable to
            specify a function to compare them
            (e.g. lambda q, s: q.split('.')[0] == s.split('|')[1])
    """
    ## parse data into dict, only keep qseqid, sseqid, and bitscore (indexed by qseqid)
    print("Parsing data")
    if fields:
        i_qseqid = fields.index("qseqid")
        i_sseqid = fields.index("sseqid")
        i_bitscore = fields.index("bitscore")
    else:
        i_qseqid = 0
        i_sseqid = 1
        i_bitscore = -1
    dat = blast6_to_qseqid_dict(rblast, i_sseqid, i_bitscore, i_qseqid = i_qseqid)
    ## get max bitscore value for each qseqid
    print("Getting max score hits")
    dat_max = {}
    for qseqid, qdat in dat.items():
        max_score = max(qdat, key = lambda x: float(x[1]))[1]
        dat_max[qseqid] = [x for x in qdat if x[1] == max_score]
    ## if targets provided, filter by that too
    if targets is not None:
        print("Filtering for targets")
        dat_filt = {}
        compared = {}
        def is_target(sseqid):
            if sseqid not in compared:
                val = False
                for target in targets:
                    if eq(sseqid, target):
                        val = True
                        break
                compared[sseqid] = val
            return compared[sseqid]
        for qseqid, qdat in dat_max.items():
            passed = 0
            for hit in qdat:
                is_a_target = is_target(hit[0])
                if is_a_target:
                    passed += 1
                    hit.append('y')
                else:
                    hit.append('n')
            if passed:
                dat_filt[qseqid] = qdat
        dat_max = dat_filt
    ## write
    print("Writing to file")
    with open(fout, "w+") as f:
        f.write('\t'.join(["qseqid", "sseqid", "bitscore"] +
                          ([] if not targets else ["sseq_is_target"])) + '\n')
        for qseqid, qdat in dat_max.items():
            for hit_dat in qdat:
                try:
                    f.write('\t'.join([qseqid] + hit_dat) + '\n')
                except Exception as e:
                    print(qseqid, hit_dat)
                    raise e
    return
