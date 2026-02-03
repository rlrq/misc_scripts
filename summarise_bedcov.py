def total_coverage(fname):
    with open(fname, 'r') as f:
        line = f.readline().rstrip().split('\t')
        curr_range = tuple(line[:3])
        curr_cov = int(line[-1])
        for line in f:
            line = line.rstrip().split('\t')
            line_range = tuple(line[:3])
            line_depth = int(line[-1])
            if line_range == curr_range:
                curr_cov += line_depth
            else:
                yield [curr_range[0], int(curr_range[1]), int(curr_range[2]), curr_cov]
                curr_range = line_range
                curr_cov = line_depth
    return

def average_depth(fname):
    for range_coverage in total_coverage(fname):
        chrom, start, end, coverage = range_coverage
        avg_depth = coverage/(end-start)
        yield [chrom, start, end, avg_depth]
    return

def write_average_depth(fin, fout):
    with open(fout, 'w+') as f:
        for entry in average_depth(fin):
            _ = f.write('\t'.join(map(str, entry)) + '\n')
    return
