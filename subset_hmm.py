#!/usr/bin/python3
import argparse

def split_to_tuple(inpt, delim = ','):
    return (tuple() if len(inpt) == 0 else tuple(inpt.split(',')))

def subset_hmm(fin, fout, names = '', accs = '', accs_v = '', exclude = ''):
    names = set(split_to_tuple(names, ','))
    accs = set(split_to_tuple(accs, ','))
    accs_v = set(split_to_tuple(accs_v, ','))
    exclude = set(map(lambda x: x.ljust(5, ' '), split_to_tuple(exclude, ',')))
    def make_keep():
        functions = []
        if names:
            functions.append(lambda name, acc: (name in names))
        if accs:
            functions.append(lambda name, acc: (acc.split('.')[0] in accs))
        if accs_v:
            functions.append(lambda name, acc: (acc in accs_v))
        return lambda name, acc: any(map(lambda f: f(name, acc), functions))
    keep = make_keep()
    entry = []
    first = True
    store = None
    name = None
    acc = None
    ## open file for reading and writing
    with open(fin, 'r') as f, open(fout, 'w+') as out:
        for line in f:
            ## if first line in entry
            if first:
                entry = [line]
                store = True
                first = False
                continue
            ## start parsing fields
            field = line[:5]
            if field == "NAME ":
                name = line[6:-1]
                entry.append(line)
            ## if at 'ACC' attribute, check if this is an entry we want to keep
            elif field == "ACC  ":
                acc = line[6:-1]
                ## if not keeping, set store=False to skip all subsequent lines
                ## (until we encounter attribute 'NAME' again)
                if not keep(name, acc):
                    store = False
                ## if keeping, store line
                else:
                    entry.append(line)
                continue
            ## store lines
            elif store:
                if field in exclude:
                    continue
                entry.append(line)
            ## end of entry, write all lines in entry if store=True
            if line == "//\n":
                if store:
                    for l in entry:
                        out.write(l)
                ## raise flag to start searching for next entry
                first = True
    return      


parser = argparse.ArgumentParser(description="subset HMM file")
parser.add_argument("--in", type=str, help="path to input HMM file for subsetting",
                    dest = "fin")
parser.add_argument("--out", type=str, help="path to output (subsetted) HMM file",
                    dest = "fout")
parser.add_argument("--name", type=str, help="comma-separated profile names to subset",
                    default = '')
parser.add_argument("--acc", type=str, help="comma-separated profile accessions to subset",
                    default = '')
parser.add_argument("--acc-v", type=str,
                    help="comma-separated profile accessions with version to subset",
                    default = '')
parser.add_argument("--exclude", type=str,
                    help="comma-separated field names to exclude",
                    default = '')
# parser.add_argument("--essential", action="store_true", default=False,
#                     help="write only essential fields")
args = parser.parse_args()

subset_hmm(args.fin, args.fout, names = args.name, accs = args.acc, accs_v = args.acc_v,
           exclude = args.exclude)
