#!/usr/bin/python3

import re

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import data_manip

filename = "/mnt/chaelab/rachelle/pav_haplotype/data/homologues/phytozome13.AT1G72930.homologs.tsv"
fout = "/mnt/chaelab/rachelle/pav_haplotype/results/homologues/phytozome13.B5.AT1G72930.homologs.tsv"
include_gid = ["AT1G17600", "AT1G17610", "AT1G17615", "AT1G17680", "AT1G72840", "AT1G72850", "AT1G72852", "AT1G72860", "AT1G72870", "AT1G72890", "AT1G72900", "AT1G72910", "AT1G72920", "AT1G72930", "AT1G72940", "AT1G72950", "AT4G09360", "AT4G09420", "AT4G09430", "AT5G40090", "AT5G40100", "AT5G48770", "AT5G48775", "AT5G48780"]
include_gid_set = set(include_gid)

tid2gid_v0 = lambda tid:re.search(".+?(?=\.\d+$)", tid).group(0)
reference_organisms_v0 = ["A.thaliana Araport11"]

def filter_clade_homologues_phytozome13(f_homologues, fout, include_gid, tid2gid = tid2gid_v0,
                                        reference_organisms = reference_organisms_v0):
    ## coerce some stuff into sets
    reference_organisms = set(reference_organisms)
    include_gid_set = set(include_gid)
    ## read & sort homologues data
    get, data = data_manip.parse_get_data(f_homologues, delim = '\t', parse_num = True)
    data_organism = [x for x in data if get(x, "Organism") in reference_organisms]
    data_organism_sorted = sorted(data_organism, key = lambda x:get(x, "Score"), reverse = True)
    ## process, starting from highest to lowest score
    include_threshold = None
    for entry in data_organism_sorted:
        if tid2gid(get(entry, "Transcript Name")) in include_gid_set:
            include_threshold = get(entry, "Score")
        else:
            break ## break upon encountering first transcript not belonging to genes to include
    ## get entries to write
    keep_entries = [x for x in data if get(x, "Score") >= include_threshold]
    ## write
    with open(fout, "w+") as f:
        _ = f.write('\t'.join(get(get_cols = True)) + '\n')
        for entry in keep_entries:
            _ = f.write('\t'.join(entry) + '\n')
    return
