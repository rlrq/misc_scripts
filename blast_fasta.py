import json
import os
import re

from Bio import Seq


def get_hits(blastj) -> list:

    f = open(blastj, 'r')
    contents = json.load(f)
    f.close()

    results = contents["BlastOutput2"]["report"]["results"]["bl2seq"]
    hits = list(map(lambda x: x["hits"], results))
    hits_out = []
    for hit in hits:
        hits_out.extend(hit)

    return hits_out


# only returns True if there is an overlap and ranges are in the same direction
def has_overlap(t1, t2):
    return ((t1[0] < t1[1] and t2[0] < t2[1]) and
            (t1[0] <= t2[0] <= t1[1] or
             t2[0] <= t1[0] <= t2[1] or
             t1[0] <= t2[1] <= t1[1] or
             t2[0] <= t1[1] <= t2[1])) or\
           ((t1[0] > t1[1] and t2[0] > t2[1]) and
            (t1[0] >= t2[0] >= t1[1] or
             t2[0] >= t1[0] >= t2[1] or
             t1[0] >= t2[1] >= t1[1] or
             t2[0] >= t1[1] >= t2[1]))


def is_ascending(t) -> bool:
    return t[0] <= t[1]


def combine_tuprange(t1, t2) -> tuple:
    if not t1:
        return t2
    elif not t2:
        return t1
    if has_overlap(t1, t2):
        if is_ascending(t1):
            return tuple((min(*t1, *t2), max(*t1, *t2)))
        else:
            return tuple((max(*t1, *t2), min(*t1, *t2)))
    raise ValueError


def combine_tupranges(*tup_ranges) -> tuple:
    tup_ranges = [x for x in tup_ranges if (x != [] and x != ())]
    output = tup_ranges[0]
    for tup_range in tup_ranges[1:]:
        output = combine_tuprange(output, tup_range)
    return output


def get_overlaps(lst, tup_range) -> list:
    output = []
    for lst_range in lst:
        if has_overlap(lst_range, tup_range):
            output.append(lst_range)
    return output


def get_hits_minlen(blastj, lower: float=0, blast2=True) -> dict:

    """
    Retrieves the ranges (combined) of hits, given as lists of tuples (where each tuple
    corresponds with one (<lower>,<higher>) subject range (inclusive) indexed by subject
    title.
    :param blastj: str. path to .json file of blast output (tblastn)
    :param lower: lowest accepted fraction of query aligned (0 <= lower <= 1)
    :return: dict. see function description.
    """

    f = open(blastj, 'r')
    contents = json.load(f)
    f.close()

    if blast2:
        results = contents["BlastOutput2"]["report"]["results"]["bl2seq"]
    else:
        results = contents["BlastOutput2"]["report"]["results"]["search"]["hits"]
    hits_out = {}

    for hits in results:

        curr_hit_ranges = []
        if blast2:
            # if no hits, move on
            if not hits["hits"]:
                continue

            # instantiate some variables
            curr_hit_name = hits["hits"][0]["description"][0]["title"]
            hit_iter = hits["hits"]
            
        else:
            hit_iter = [hits]
            curr_hit_name = hits["description"][0]["title"]

        # iterate through hits
        for hits_entries in hit_iter:
            for hit in hits_entries["hsps"]:
                if blast2:
                    ref_len = hits["query_len"]
                else:
                    ref_len = hits["len"]
                if hit["align_len"]/ref_len >= lower:

                    # get subject range in alignment
                    curr_range = tuple(sorted((hit["hit_from"], hit["hit_to"])))
                    if not curr_hit_ranges:
                        curr_hit_ranges.append(curr_range)
                        continue

                    # find any overlaps and combine
                    overlaps = get_overlaps(curr_hit_ranges, curr_range)
                    for overlap in overlaps:
                        curr_hit_ranges.remove(overlap)
                    curr_hit_ranges.append(combine_tupranges(curr_range, *overlaps))
                    
        if curr_hit_ranges:
            hits_out[curr_hit_name] = curr_hit_ranges

    return hits_out


def get_hits_minlen_diraware(blastj, lower: float=0) -> dict:

    """
    Retrieves the ranges (combined) of hits, given as lists of tuples (where each tuple
    corresponds with one (<lower>,<higher>) subject range (inclusive) indexed by subject
    title.
    :param blastj: str. path to .json file of blast output (tblastn)
    :param lower: lowest accepted fraction of query aligned (0 <= lower <= 1)
    :return: dict. see function description.
    """

    f = open(blastj, 'r')
    contents = json.load(f)
    f.close()

    results = contents["BlastOutput2"]["report"]["results"]["bl2seq"]
    hits_out = {}

    for hits in results:

        # if no hits, move on
        if not hits["hits"]:
            continue

        # instantiate some variables
        curr_hit_ranges = []
        curr_hit_name = hits["hits"][0]["description"][0]["title"]

        # iterate through hits
        for hits_entries in hits["hits"]:
            for hit in hits_entries["hsps"]:
                if hit["align_len"]/hits["query_len"] >= lower:

                    # get subject range in alignment
                    curr_range = tuple((hit["hit_from"], hit["hit_to"]))

                    plus_strand = bool(hit["hit_frame"] > 0)

                    if not plus_strand:
                        curr_range = sorted(curr_range, reverse=True)

                    if not curr_hit_ranges:
                        curr_hit_ranges.append(curr_range)
                        continue

                    # find any overlaps and combine
                    overlaps = get_overlaps(curr_hit_ranges, curr_range)
                    for overlap in overlaps:
                        curr_hit_ranges.remove(overlap)
                    curr_hit_ranges.append(combine_tupranges(curr_range, *overlaps))

        if curr_hit_ranges:
            hits_out[curr_hit_name] = curr_hit_ranges

    return hits_out


def get_fastaseq_from_range(fname, ranges, fout=None, linelen=60, translate=False):

    # read fasta file
    fasta_seqs = read_fasta(fname)
    output = {}

    # iterate through tuple ranges in range and write to output dict
    # iterate through all entries in fasta file
    for name, seq in fasta_seqs.items():

        # ignore entry if not found in key of ranges
        for tup_range in ranges.get(name, ()):

            # obtain base name
            name_base = (name if "NODE" not in name else
                         re.findall(re.compile("(.+?)_cov_.+"), name)[0])
            if is_ascending(tup_range):
                output_seq = seq[tup_range[0]-1:tup_range[1]]
            else:
                output_seq = Seq.reverse_complement(seq[tup_range[1]-1:tup_range[0]])
            if translate:
                output_seq = Seq.translate(output_seq)

            output["{}_{}-{}".format(name_base, *tup_range)] = output_seq

    # write output file
    dict_to_fasta(output,
                  (fout if fout else os.path.splitext(fname)[0]+"_range.fasta"),
                  linelen=linelen)


def translate_fasta(fname, fout=None, linelen=60):

    f = read_fasta(fname)
    for name, seq in f.items():
        f[name] = Seq.translate(seq)
    dict_to_fasta(f,
                  fout if fout else os.path.splitext(fname)[0]+"_translated.fasta",
                  linelen=linelen)


def get_hit_titles(hits) -> list:

    return list(set(map(lambda x: x["description"][0]["title"], hits)))


def read_fasta(fname) -> dict:

    """
    Corrals sequences in fasta file into {<seqname>: <sequence>} dictionary
    :param fname: str. file name
    :return: dictionary
    """

    f = open(fname, 'r')
    seqs = {}
    curr_seqname = ''
    curr_seq = ''
    for line in f.readlines():
        if line[0] == '>':
            if curr_seqname and curr_seq:
                seqs[curr_seqname] = curr_seq
            curr_seqname = line[1:-1]
            curr_seq = ''
        else:
            curr_seq += line[:-1]
    seqs[curr_seqname] = curr_seq
    f.close()
    return seqs


def break_fasta(fname, n, linelen=60):

    f = read_fasta(fname)
    name, ext = os.path.splitext(fname)
    seqs = list(f.items())
    for i in range(0, len(seqs), n):
        tup_to_fasta(seqs[i:i+n], name + '_' + str((i/n) + 1) + ext, linelen=linelen)


def tup_to_fasta(t, fout, linelen=60):

    f = open(fout, 'w+')
    output = ''
    for name, seq in t:
        output += '>' + name + '\n'
        output += '\n'.join([seq[i:i + linelen] for i in range(0, len(seq), linelen)]) + '\n'
    f.write(output)
    f.close()


def dict_to_fasta(d, fout, linelen=60):

    f = open(fout, 'w+')
    output = ''
    for name, seq in d.items():
        output += '>' + name + '\n'
        output += '\n'.join([seq[i:i+linelen] for i in range(0, len(seq), linelen)]) + '\n'
    f.write(output)
    f.close()


def filter_fasta(fname, seq_names, outname=None, linelen=60) -> None:

    faf = read_fasta(fname)
    output = ''
    for name, seq in faf.items():
        if name in seq_names:
            output += '>' + name + '\n'
            output += '\n'.join([seq[i:i+linelen] for i in range(0, len(seq), linelen)]) + '\n'
    outname = (os.path.splitext(outname)[0] if outname else (fname + "_filtered")) + ".fasta"
    outf = open(outname, 'w+')
    outf.write(output)
    outf.close()

    return


def filter_blast_fasta(blast_file, fa_file, outname=None, linelen=60, gethits=get_hits):

    hit_names = get_hit_titles(gethits(blast_file))
    filter_fasta(fa_file, hit_names, outname=outname, linelen=linelen)


def filter_blast_fasta_multi(blast_files, fa_file, linelen=60, gethits=get_hits):

    fa_name, fa_ext = os.path.splitext(fa_file)
    for blast_file in blast_files:
        blast_name, blast_ext = os.path.splitext(blast_file)
        filter_blast_fasta(blast_file, fa_file,
                           outname=fa_name+'_'+blast_name+"_filtered.fasta",
                           linelen=linelen, gethits=gethits)


blast_files_ey152 = ("tblastn_RPP5_NBS_EY15-2.json",
                     "tblastn_RPP5_TIR_EY15-2.json",
                     "tblastn_SNC1_NBS_EY15-2.json",
                     "tblastn_SNC1_TIR_EY15-2.json")

blast_files_ice228 = ("tblastn_RPP5_NBS_ICE228.json",
                      "tblastn_RPP5_TIR_ICE228.json",
                      "tblastn_SNC1_NBS_ICE228.json",
                      "tblastn_SNC1_TIR_ICE228.json")

snc1_nbs = "pseaamieelaedvlrktmtpsddfgdlvgienhieaiksvlcleskearimvgiwgqsgigkstigralysklsiqfhhrafitykstsgsdvsgmklrwekellseilgqkdikiehfgvveqrlkqqkvlillddvdsleflktlvgkaewfgsgsriivitqdrqllkaheidliyevefpsehlaltmlcrsafgkdsppddfkelafevaklagnlplglsvlgsslkgrtkewwmemmprlrnglngdimktlrvsydrlhqkdqdmflyiaclfngfevsyvkdllkdnvgftmlteksliritpdgyiemhnlleklgreidrakskgnpgkrrfltnfedihevvtekt"
snc1_tir = "mmdtskdddmeiasssgsrrydvfpsfrgedvrdsflshllkelrgkaitfiddeiersrsigpellsaikesriaivifsknyasstwclnelveihkcytnlnqmvipiffhvdasevkkqtgefgkvfeetckaksedekqswkqalaavavmagydlrk"
rpp5_nbs = "pneahmvekisndvsnklitrskcfddfvgieahieaiksvlcleskearmvgiwgqsgigkstigralfsqlsiqfplrafltykstsgsdvsgmklswekellseilgqkdikiehfgvveqrlkhkkvlillddvdnlefkltlvgkaewfgsgsriivitqdrqflkahdidlvyevklpsqglaltmlcrsafgkdsppddfkelafevaklaghlplglnvlgsslrrrgkkewmemmprlrnglngdimktlrvsydrlhqkdqdmflciaclfngfevsyvkdllednvgltmlseksliritpdghiemhnlleklgreidrakskgnpgkrqfltnfedihevvtekt"
rpp5_tir = "maassssgrrrydvfpsfsgvdvrktflshliealdgksintfidhgiersrtiapelisairearisivifsknyasstwclnelveihkcfndlgqmvipvfydvdpsevrkqtgefgkvfektcevskdkqpgdqkqrwvqaltdianiagedlln"
