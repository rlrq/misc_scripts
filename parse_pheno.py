import sys
sys.path.append("/home/rachelle/scripts")

from basics import *
from data_manip import *

# f_acc_map = "/home/rachelle/data/1001genomes/1135acc.csv"
# get_acc, dat_acc = parse_get_data(f_acc_map)
# acc_map = {str(get_acc(x, "tg_ecotypeid")): get_acc(x, "CS_number") for x in dat_acc}

def make_acc_map():
    get_acc, dat_acc = parse_get_data("/mnt/chaelab/rachelle/data/1001genomes/1135acc.csv")
    acc_map = {str(get_acc(x, "tg_ecotypeid")): get_acc(x, "CS_number") for x in dat_acc}
    return acc_map

## parse hap file to dict
def hap_to_dict(fname, acc_map = make_acc_map(), valid_accs = None):
    get, dat_hap = parse_get_data(fname)
    valid_accIDs = sorted(set(x for x in (','.join([str(y) for y in get(dat_hap, "accIDs")])).split(',')
                              if (acc_map is None or x in acc_map)))
    haps = {accID: [] for accID in valid_accIDs}
    hap_names = sorted(set(get(dat_hap, "domain")))
    ## convert data to 0 and 1
    for entry in sorted(dat_hap):
        accIDs = str(get(entry, "accIDs")).split(',')
        to_append = 1 if get(entry, "group") == 1 else 0
        for accID in accIDs:
            if accID in haps:
                haps[accID].append(to_append)
    return haps, hap_names

def clearHap_to_dict(fname, acc_map = make_acc_map()):
    get, dat_hap = parse_get_data(fname)
    ## filter accs
    valid_accIDs = sorted(set(x for x in (','.join([str(y) for y in get(dat_hap, "accIDs")])).split(',')
                              if (acc_map is None or x in acc_map)))
    ## parse dat
    haps = {accID: [] for accID in valid_accIDs}
    hap_names = []
    ## convert data to 0 and 1
    for entry in sorted(dat_hap):
        if get(entry, "group") == "missing": continue
        hap_name = '.'.join(str(x) for x in get(entry, "domain", "group"))
        hap_names.append(hap_name)
        accIDs = str(get(entry, "accIDs")).split(',')
        for accID in haps:
            if accID in accIDs: haps[accID].append(1)
            else: haps[accID].append(0)
    return haps, hap_names

## write hap dictionary to file
def write_plink_hap(fout, haps, hap_names):
    with open(fout, "w+") as f:
        header = ["FID", "IID"] + hap_names
        data = [[k, k] + v for k, v in haps.items()]
        f.write('\n'.join(['\t'.join(map(str, x)) for x in [header] + data]))
    return

## see f_hap for input file format
def plink_from_hap(fin, fout, acc_map = make_acc_map(), get = None):
    write_plink_hap(fout, *hap_to_dict(fin, acc_map = acc_map))
    return

## see f_hap for input file format
def plink_from_clearhap(fin, fout, acc_map = make_acc_map(), valid_accs = None):
    write_plink_hap(fout, *clearHap_to_dict(fin, acc_map = acc_map, valid_accs = valid_accs))
    return

## takes in dictionary of filenames containing accIDs ('\n'-separated) to be assigned 1 (binary value),
## # where vdw accIDs not in the list are assigned 0
def plink_from_accs(d_fhaps, all_acc, acc_map, fout):
    d_acc = {k: splitlines(v) for k, v in d_fhaps.items()}
    valid_accIDs = set(x for x in all_acc if x in acc_map)
    haps = {accID: [] for accID in sorted(valid_accIDs)}
    hap_names = sorted(d_acc.keys())
    ## convert data to 0 and 1
    for hap_name in hap_names:
        accIDs = set(d_acc[hap_name])
        for accID, hap in haps.items():
            hap.append(1 if accID in accIDs else 0)
    ## write
    write_plink_hap(fout, haps, hap_names)
    return

