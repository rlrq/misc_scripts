
## executes blast w/ outfmt 6 and prepends colnames to file
## header should either be an ordered iterable of fields or str where fields are comma-separated
def blast6(blastf, header, fout, **kwargs):
    ## parse header
    if isinstance(header, str):
        header = header.split(',')
    ## execute blast
    blast_cline = blastf(out = fout, outfmt = f"'6 {' '.join(header)}'", **kwargs)
    blast_cline()
    ## reformat output file
    with open(fout, 'r') as f:
        dat = f.read()
    with open(fout, "w+") as f: ## prepend header
        f.write('\t'.join(header) + '\n' + dat)
    return
