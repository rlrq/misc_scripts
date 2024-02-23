#!/usr/bin/python3

def phylip_to_long(fname, fout, transform = lambda x:x, filter = lambda x:True):
    with open(fname, 'r') as fin:
        _ = fin.readline()
        names = [line.split()[0] for line in fin]
    with open(fname, 'r') as fin:
        _ = fin.readline()
        with open(fout, "w+") as f:
            for line in fin:
                vals = line.split()
                name = vals[0]
                vals = [transform(x) for x in vals[1:]]
                for i, val in enumerate(vals):
                    if filter(val):
                        f.write(f"{name}\t{names[i]}\t{val}\n")
    return
            
# def wide_to_long_fastprot(fname, fout, transform = lambda x:x, filter = lambda x:True):

def phylip_to_long_keep_min(fname, fout, transform = lambda x:x, filter = lambda x:True):
    with open(fname, 'r') as fin:
        _ = fin.readline()
        mat = [line.split() for line in fin]
        names = [row[0] for row in mat]
        mat = [[transform(val) for val in row[1:]] for row in mat]
    with open(fout, "w+") as f:
        for i in range(len(names)-1):
            for j in range(i, len(names)):
                try:
                    val = min(mat[i][j], mat[j][i])
                except Exception as e:
                    print(i, j, names[i], names[j])
                    raise e
                if filter(val):
                    f.write(f"{names[i]}\t{names[j]}\t{val}\n")
    return

def long_keep_min(fname, fout, transform = lambda x:float(x), filter = lambda x:True):
    ## get matrix size
    with open(fname, 'r') as fin:
        names = list(set(line.split()[0] for line in fin))
    ## generate name-index map
    names_i = {name: i for i, name in enumerate(names)}
    ## generate matrix
    mat = [[None]*len(names) for i in names]
    ## store transformed values in matrix
    with open(fname, 'r') as fin:
        for line in fin:
            nameA, nameB, val = line.split()
            mat[names_i[nameA]][names_i[nameB]] = transform(val)
    ## write only minimum
    with open(fout, "w+") as f:
        for i in range(len(names)-1):
            for j in range(i+1, len(names)):
                val = min(mat[i][j], mat[j][i])
                if filter(val):
                    f.write(f"{names[i]}\n{names[j]}\t{val}\n")
    return
