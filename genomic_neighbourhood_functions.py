#!/usr/bin/python3

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import basics
import data_manip as dm
import annotation as ann

## outputs a dictionary {'<chrom>': [<bed entries as list, sorted by start then end positions>]}
## fields for first start and end positions (index 1 and 2) will be converted to integers
def read_bed_v1(fname):
    entries = [x.split('\t') for x in dm.splitlines(fname)]
    ## group entries by chrom
    output = {}
    for entry in entries:
        ## add entry; convert start and end positions to integers
        chrom = entry[0]
        if chrom not in output:
            output[chrom] = []
        output[chrom].append([entry[0], int(entry[1]), int(entry[2])] + entry[3:])
    ## sort output
    output = {chrom: sorted(entries, key = lambda e:e[1:3])
              for chrom, entries in output.items()}
    return output

## read gff file
## outputs a dictionary {'<chrom>': [<bed entries as list, sorted by start then end positions>]}
## fields for first start and end positions (index 1 and 2) will be converted to integers
def read_gff_v1(fname):
    gff = ann.GFF(fname = fname, fmt = "GFF", memsave = True)
    ## group entries by chrom
    output = {}
    for entry in gff:
        chrom = entry.molecule
        if chrom not in output:
            output[chrom] = []
        output[chrom].append(entry)
    ## sort output
    sort_by_start = {chrom: sorted(entries, key = lambda e:(e.start,e.end))
                     for chrom, entries in output.items()}
    sort_by_end = {chrom: sorted(entries, key = lambda e:(e.end,e.start))
                   for chrom, entries in output.items()}
    return sort_by_start, sort_by_end

def format_overlap_v1(features):
    output = []
    for entry in features:
        ids = entry.get_attr("ID", fmt = list)
        ids = ids if len(ids) != 0 else ['.']
        parents = entry.get_attr("Parent", fmt = list)
        parents = parents if len(parents) != 0 else ['.']
        for _id in ids:
            for parent in parents:
                output.append({"chrom": entry.molecule, "start": entry.start, "end": entry.end,
                               "feature_type": entry.type, "ID": _id, "Parent": parent,
                               "distance": '.',
                               "position": "overlap", "proximity_rank": '.'})
    return output

def format_upstream_v1(features, start_pos):
    output = []
    for i, entry in enumerate(features[-1::-1]):
        ids = entry.get_attr("ID", fmt = list)
        ids = ids if len(ids) != 0 else ['.']
        parents = entry.get_attr("Parent", fmt = list)
        parents = parents if len(parents) != 0 else ['.']
        for _id in ids:
            for parent in parents:
                output.append({"chrom": entry.molecule, "start": entry.start, "end": entry.end,
                               "feature_type": entry.type, "ID": _id, "Parent": parent,
                               "distance": (start_pos - entry.end + 1),
                               "position": "upstream", "proximity_rank": i+1})
    return output

def format_downstream_v1(features, end_pos):
    output = []
    for i, entry in enumerate(features):
        ids = entry.get_attr("ID", fmt = list)
        ids = ids if len(ids) != 0 else ['.']
        parents = entry.get_attr("Parent", fmt = list)
        parents = parents if len(parents) != 0 else ['.']
        for _id in ids:
            for parent in parents:
                output.append({"chrom": entry.molecule, "start": entry.start, "end": entry.end,
                               "feature_type": entry.type, "ID": _id, "Parent": parent,
                               "distance": (entry.start - end_pos + 1),
                               "position": "downstream", "proximity_rank": i+1})
    return output

def get_neighbourhood_v1(start, end, gff_by_start, gff_by_end,
                         feature_type = ["gene", "pseudogene"], num_flank = 1, max_distance = 10000):
    feature_type = set(feature_type)
    overlap = [e for e in gff_by_start
               if basics.has_overlap((start, end), (e.start, e.end), end_inclusive = False)]
    upstream = [e for e in gff_by_end if start >= e.end and e.type in feature_type]
    downstream = [e for e in gff_by_start if end <= e.start and e.type in feature_type]
    upstream_subset = upstream[-num_flank:]
    downstream_subset = downstream[:num_flank]
    return [e for e in (format_upstream_v1(upstream_subset, start) +
                        format_overlap_v1(overlap) +
                        format_downstream_v1(downstream_subset, end))
            if (e["distance"] == '.' or e["distance"] <= max_distance)]

def get_neighbourhoods_v1(grouped_bed, gff_by_start, gff_by_end, **kwargs):
    output = []
    for chrom, entries in grouped_bed.items():
        gff_by_start_subset = gff_by_start.get(chrom, [])
        gff_by_end_subset = gff_by_end.get(chrom, [])
        for entry in entries:
            output.append({"bed_entry": entry,
                           "neighbourhood": get_neighbourhood_v1(entry[1], entry[2],
                                                                 gff_by_start_subset, gff_by_end_subset,
                                                                 **kwargs)})
    return output

default_neighbour_fields_v1 = ["position", "proximity_rank", "distance",
                               "chrom", "start", "end", "feature_type", "ID", "Parent"]
def write_neighbours_v1(fout, neighbourhoods, bed_columns = None, index = 1, end_inclusive = True,
                        neighbour_fields = default_neighbour_fields_v1):
    if bed_columns is None:
        bed_columns = list(range(max(len(e["bed_entry"]) for e in neighbourhoods)))
    with open(fout, 'w+') as f:
        header = neighbour_fields + [f"bed_{i+1}" for i in range(len(bed_columns))]
        _ = f.write('\t'.join(header) + '\n')
        for entry in neighbourhoods:
            bed_entry = entry["bed_entry"]
            neighbourhood = entry["neighbourhood"]
            bed_to_write = [bed_entry[i] for i in bed_columns]
            for neighbour in neighbourhood:
                try:
                    neighbour["start"] = neighbour["start"] + index
                    neighbour["end"] = neighbour["end"] + index - (1 if end_inclusive else 0)
                except Exception as e:
                    print(neighbour)
                    raise e
                to_write = [neighbour.get(field, '.') for field in neighbour_fields] + bed_to_write
                _ = f.write('\t'.join(map(str, to_write)) + '\n')
    return

def get_neighbourhoods_from_files_v1(f_bed, f_gff, f_out,
                                     feature_type = ["gene", "pseudogene"], num_flank = 1, max_distance = 10000):
    grouped_bed = read_bed_v1(f_bed)
    gff_by_start, gff_by_end = read_gff_v1(f_gff)
    neighbourhoods = get_neighbourhoods_v1(grouped_bed, gff_by_start, gff_by_end,
                                           feature_type = feature_type, num_flank = num_flank,
                                           max_distance = max_distance)
    write_neighbours_v1(f_out, neighbourhoods, index = 1, end_inclusive = True,
                        neighbour_fields = default_neighbour_fields_v1)
    return
