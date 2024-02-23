def make_custom_get(header, parse_num = True):
    def get_col(colname, data):
        return [get_col_in_row(x, colname) for x in data]
    def get_col_in_row(row, colname):
        if not colname in header:
            print(f"Column '{colname}' is not found in headers {header}")
        output = row[header.index(colname)]
        if isinstance(output, (list, tuple)):
            if [str(x).isdigit() for x in output].count(True) == len(output):
                output = [int(x) for x in output]
            elif [str(x).replace('.','',1).replace('-','',1).isdigit() \
                  for x in output].count(True) == len(output):
                output = [float(x) for x in output]
            return output
        else:
            return output if (parse_num is False
                              or (not str(output).replace('.','',1).replace('-','',1).isdigit())) else \
                float(output) if not str(output).isdigit() else int(output)
    def helper(data = None, *colnames, get_cols = False, suppress_print = False, ncol = False,
               return_list = False):
        if get_cols:
            return header
        if ncol:
            return len(header)
        if len(data) == 0:
            if not suppress_print:
                print("No data found; returning empty list")
            return []
        if isinstance(data[0], (list, tuple)):
            output = [get_col(colname, data) for colname in colnames]
            return output[0] if len(output) == 1 else [[output[r][c] for r in range(len(output))] \
                                                       for c in range(len(output[0]))]
        else:
            output = [get_col_in_row(data, colname) for colname in colnames]
            return output[0] if (len(output) == 1 and not return_list) else output
    return helper

def parse_get_data(fname, delim = '\t', detect = True, parse_num = True):
    if detect:
        ext = fname.split('.')[-1]
        if ext == "csv":
            delim = ','
        elif ext == "tsv":
            delim = '\t'
    # data = [(''.join([c for c in line if c != '\n'])).split(delim)
    #         for line in open(fname, 'r').readlines()]
    data = [line.split(delim) for line in splitlines(fname)]
    get = make_custom_get(data[0], parse_num = parse_num)
    return get, data[1:]

def extend_get(get_f, extension):
    if callable(extension):
        return make_custom_get(get_f(None, get_cols = True) +
                               extension(None, get_cols = True))
    else:
        return make_custom_get(type(extension)(get_f(None, get_cols = True)) +
                               extension)

def splitlines(fname, ignore_empty_lines = True):
    data = open(fname, 'r').read().split('\n')
    if ignore_empty_lines:
        return [line for line in data if line]
    else:
        return data

def filter_dict(d, keys):
    keys = set(keys)
    return {k: v for k, v in d.items() if k in keys}

def write_delim(fname, dat, delim = '\t', coerce = True):
    open(fname, "w+").write('\n'.join([delim.join([str(x) if coerce else x for x in line])
                                       for line in dat]) + '\n')
    return

def sort_get_multi(iterable, get, *cols): ## the order is important!
    ## assumes conserved order sort
    for col in cols:
        iterable.sort(key = lambda x: get(x, col))
    return
