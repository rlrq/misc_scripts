## variables format: [(<varname>, <value>, <shortname>, [<longname1>, <longname2>, ...])]
## ignore format: [<varname1>, <varname2>, ...]
def write_log(fout, variables = [], header = '', ignore = [], append = False):
    output = []
    for var in [var for var in variables if var[0] not in ignore]:
        output.append(['|'.join(list(var[2]) + list(var[3])), str(var[1])])
    to_write = '\n'.join(['\t'.join(entry) for entry in output]) + '\n'
    if header:
        to_write = header + '\n' + to_write
    if append:
        open(fout, "a+").write(to_write)
    else:
        open(fout, "w+").write(to_write)
    return
