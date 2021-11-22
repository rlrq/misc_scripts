def index_file(fname, chunk_lines = 10000):
    output = [0]
    with open(fname, 'r') as f:
        i = 0
        while True:
            i += 1
            line = f.readline()
            if not line:
                break
            if i % chunk_lines == 0:
                output.append(f.tell())
    return output

def get_line(i, fname, index, chunk_lines):
    i_chunk = i // chunk_lines
    mod_chunk = i % chunk_lines
    with open(fname, 'r') as f:
        f.seek(index[i_chunk])
        for i, line in enumerate(f):
            if i == mod_chunk: break
    return line

class IndexedFile:
    def __init__(self, filename, chunk_lines = 10000, skip = lambda x: False):
        self.filename = filename
        self.chunk_lines = chunk_lines
        self.indices = None
        ## if skip(line) returns True, the line will not be indexed
        ##  and will be ignored (and will not add to line count) when using get_line
        self.skip = skip
        self.index()
    def __iter__(self):
        with open(self.filename, 'r') as f:
            for line in f:
                if self.skip(line): continue
                else: yield
    def index(self, chunk_lines = None):
        if chunk_lines is not None:
            self.chunk_lines = chunk_lines
        indices = [0]
        with open(self.filename, 'r') as f:
            i = 0
            while True:
                line = f.readline()
                if not line: break
                if self.skip(line): continue
                i += 1
                if i % self.chunk_lines == 0:
                    indices.append(f.tell())
        self.indices = indices
    def get_line(self, *indices, strip_newline = False):
        single_i = len(indices) == 1
        ## remove invalid indices (<0)
        is_valid = lambda i: i >= 0
        invalid = [i for i in indices if not is_valid(i)]
        if invalid: print(f"Invalid indices (will be ignored): {','.join(map(str, sorted(invalid)))}")
        indices = [i for i in indices if is_valid(i)]
        ## sort indices and split into bins according to self._chunk_lines
        bins = {}
        for i in sorted(indices):
            i_chunk = i // self.chunk_lines
            mod_chunk = i % self.chunk_lines
            bins[i_chunk] = bins.get(i_chunk, []) + [mod_chunk]
        ## jump to each bin and iterate within it
        output = []
        with open(self.filename, 'r') as f:
            for i_chunk, mod_chunks in bins.items():
                last_mod_chunk = mod_chunks[-1]
                mod_chunks = set(mod_chunks)
                ## jump to first line of bin
                f.seek(self.indices[i_chunk])
                i = 0
                for line in f:
                    if self.skip(line): continue
                    if i in mod_chunks:
                        output.append(line)
                        ## exit bin if all desired lines in bin have been visited
                        if i == last_mod_chunk: break
                    i += 1
        if strip_newline: output = [line.replace('\n', '') for line in output]
        return output if not single_i else output[0]
