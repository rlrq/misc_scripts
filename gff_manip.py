import re
import csv
import sys
import itertools

sys.path.append("/mnt/chaelab/rachelle/src")
from basics import get_recursively
from data_manip import splitlines
from index_file import IndexedFile

## parse_sep_dict and parse_attr_mod_sdict are copied from MINORg's parse_config.py
def parse_sep_sdict(s: str, item_sep: str = ';', kv_sep = ':'):
    '''
    Converts raw '<key1><kv_sep><v1>,<v2>,...,<v5><item_sep><key2><kv_sep><v6>,...,<v9>' string into 
    {'<key1>':'<v1>,<v2>,...,<v5>', '<key2>':'<v6>,...,<v9>'} dict.
    '''
    return dict( ( item.split(kv_sep)[0], kv_sep.join(item.split(kv_sep)[1:]) )
                 for item in s.split(item_sep) if item )

def parse_attr_mod_sdict(s: str, attr_sep: str = ',', feature_sep: str = ';', fa_sep: str = ':',
                         alias_sep: str = '='):
    """
    Parse attribute modifications from string to dictionary.

    Arguments:
        s (str): required, string of attribute modifications in format
            '<feature type>:<standard attribute field name>=<nonstandard attribute field name>,<standard attribute field name>=<nonstandard attribute field name>;<feature type>:<standard attribute field name>=<nonstandard attribute field name>'
            (e.g. 'mRNA:Parent=Locus_id')
        attr_sep (str): delimiter for attribute modifications of same feature type (default=',')
        feature_sep (str): delimiter for feature types (default=';')
        fa_sep (str): delimiter between feature type and attribute modifications (default=':')
        alias_sep (str): delimiter between standard attribute field name and non-standard attribute field name
            (default='=')

    Returns
    -------
    dict
        parsed attribute modifications in format
            {<feature>: {<standard attribute field name>: <nonstandard attribute field name>}}
    """
    if not s: return {}
    parse_attr_mod = lambda mappings: dict(mapping.split(alias_sep) for mapping in mappings.split(attr_sep))
    feature_dict = parse_sep_sdict(s, kv_sep = fa_sep, item_sep = feature_sep)
    return {feature: parse_attr_mod(mappings) for feature, mappings in feature_dict.items()}

def append_geneID(gff_bed_path, fout, names_col=-1, get_unique=True):

    """
    :param gff_bed_path: path to bed formatted gff file
    :names_col: column from which to extract gene ID
    """
    gff = [] # list to store all gffs

    with open(gff_bed_path, "r") as f:
        lines = f.readlines()
        print(len(lines))
        for i,v in enumerate(lines):
            temp = v[:-1].split("\t")
            geneID = re.findall("=(AT.+?G\d+)",temp[names_col])
            # if no gene ID found, skip this line
            if geneID:
                # use if replacing names column with gene ID
                temp[-1] = geneID[0]
                # use if appending gene ID to row
                # temp.append(geneID[0])
                gff.append('\t'.join(temp))

    # write features to file
    import os
    f = open(fout, "w+")
    f.write('\n'.join(gff))
    f.close()
    # sort for thoroughness
    if get_unique:
        os.system("sort -u {} | bedtools sort -i - > temp.txt; mv temp.txt {}".format(fout, fout))
    else:
        os.system("bedtools sort -i {} > temp.txt; mv temp.txt {}".format(fout, fout))

    return 0

def filter_geneIDs(fname, fout, fgeneIDs, id_sep=',', gene_col=None):

    ## if gff doesn't have geneID appended, append it
    if gene_col == None:
        append_geneID(fname, fout, gene_col=-1)
        fname = fout

    ## get list of geneIDs to be selected
    geneIDs = []
    with open(fgeneIDs) as f:
        lines = [x[:-1] for x in f.readlines()]
        for line in lines:
            geneIDs.extend(l.split(id_sep))

    gff_sel = []
    ## start filtering
    with open(fout) as f:
        lines = [x[:-1] for x in f.readlines()]
        gff_sel = ['\t'.join(x) for x in lines if (x[gene_col] in geneIDs)]

    # write features to file
    import os
    f = open(fout, "w+")
    f.write('\n'.join(gff_sel))
    f.close()
    # sort for thoroughness
    if get_unique:
        os.system("sort -u {} | bedtools sort -i - > temp.txt; mv temp.txt {}".format(fout, fout))
    else:
        os.system("bedtools sort -i {} > temp.txt; mv temp.txt {}".format(fout, fout))

    return 0

class GFF:
    
    def __init__(self, fname = None, data = [], attr_mod = {}, fmt = None, quiet = False,
                 memsave = False, chunk_lines = 1000, **kwargs):
        self._fmt = (fmt if fmt
                     else None if not fname
                     else ("BED" if fname.split('.')[-1].upper() == "BED" else "GFF3"))
        self._fname = fname
        self._data = data
        self._attr_mod = attr_mod if isinstance(attr_mod, dict) else parse_attr_mod_sdict(attr_mod)
        self._attr_fields = {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                                     "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                                     "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                                     "Is_circular": "Is_circular"}}
        self._attr_fields_inv = {}
        self._quiet = quiet
        self._kwargs = kwargs
        self._molecules = {}
        self._chunk_lines = chunk_lines
        self._indexed_file = None
        self.update_attr_fields()
        if not memsave:
            self.parse()
        elif self.has_file():
            self.index()
    
    def __iter__(self):
        if self._data:
            for entry in self._data:
                yield entry
        else:
            mk_annotation = self.make_annotation_from_str_gen(strip_newline = True)
            last_entry = None
            for entry in self.iter_raw():
                yield mk_annotation(entry)
        return
    
    def __len__(self):
        return len(self._data)
    
    def is_bed(self):
        return self._fmt.upper() == "BED"
    def is_gff(self):
        return self._fmt.upper() in {"GFF", "GFF3"}
    def has_file(self):
        return self._fname is not None
    def index(self, chunk_lines = None):
        if chunk_lines:
            self._chunk_lines = chunk_lines
        if self.has_file():
            self._indexed_file = IndexedFile(self._fname, chunk_lines = self._chunk_lines,
                                             skip = lambda line: (tuple(line[:1]) == ('#',) or line == '\n'))
    def is_indexed(self):
        return bool(self._indexed_file)
    
    def make_annotation_from_str_gen(self, strip_newline = True):
        if strip_newline:
            split_line = lambda entry: entry.replace('\n', '').split('\t')
        else:
            split_line = lambda entry: entry.split('\t')
        if self.is_bed():
            def mk_annotation(entry):
                entry = split_line(entry)
                gff_fmt = [entry[0], entry[6], entry[7], str(int(entry[1]) + 1), entry[2],
                           entry[4], entry[5], entry[8], entry[9]]
                return Annotation(gff_fmt, self, **self._kwargs)
        else: ## GFF3
            def mk_annotation(entry):
                entry = split_line(entry)
                return Annotation(entry, self, **self._kwargs)
        return mk_annotation
    
    def iter_raw(self):
        if self._indexed_file is not None:
            for entry in self._indexed_file:
                yield entry
        elif self._fname is not None:
            with open(self._fname, 'r') as f:
                for entry in f:
                    if tuple(entry[:1]) == ('#',) or entry == '\n': continue
                    yield entry
        return
    
    def parse(self):
        if self._fname is not None:
            ## start parsing
            self._data = []
            for entry in self:
                self.add_entry(entry)
            self.index_molecules()
    
    def read_file(self, fname):
        self._fname = fname
        self.parse()
    
    def index_molecules(self, reset = True):
        if reset:
            self._molecules = {}
        index = 0
        curr = None
        for entry in self:
            if entry.molecule != curr and entry.molecule not in self._molecules:
                curr = entry.molecule
                index += 1
                self._molecules[entry.molecule] = index
        return
    
    def sort(self):
        self._data.sort(key = lambda entry: (self._molecules[entry.molecule], entry.start, -entry.end,
                                             entry.source, entry.feature, entry.attributes))
    
    def add_entry(self, gff_entry):
        self._data.append(gff_entry)
    
    def update_attr_fields(self, attr_mod = None):
        '''
        Update attribute fields w/ user-provided attribute field modification dictionary.
        (Otherwise, use self._attr_fields)
        '''
        if attr_mod is None: attr_mod = self._attr_mod
        for feature, mods in attr_mod.items():
            if feature in self._attr_fields:
                for canon, atypical in mods.items():
                    self._attr_fields[feature][canon] = atypical
            else:
                self._attr_fields[feature] = mods
        self._attr_fields_inv = self.invert_attr_fields()
        return
    
    def invert_attr_fields(self):
        output = {feature: {v: k for k, v in attr_map.items()}
                 for feature, attr_map in self._attr_fields.items()}
        return output
    
    def get_subfeatures(self, feature_ids, *features, index = False):
        '''
        Get subfeatures of features w/ user-provided feature_ids.
        '''
        features = set(features)
        if type(feature_ids) is str:
            feature_ids = {feature_ids}
        else:
            feature_ids = set(feature_ids)
        indices = [i for i, entry in enumerate(self) if
                   ((not features or entry.feature in features) and
                    entry.has_attr("Parent", feature_ids))]
        if index: return indices
        else: return self.get_i(indices)
    
    def get_subfeatures_full(self, feature_ids, *features, index = False):
        '''
        Get all features that are subfeatures of user-provided feature_ids
        AND subfeatures of those subfeatures, until there are no sub-sub...sub-features left.
        '''
        features = set(features)
        output_indices = []
        curr_indices = []
        iteration_n = 0
        while True:
            iteration_n += 1
            print(f"Executing iteration {iteration_n}")
            curr_indices = self.get_subfeatures(feature_ids, index = True)
            feature_ids = set(itertools.chain(*[self.get_i(i).get_attr("ID") for i in curr_indices]))
            output_indices.extend(curr_indices)
            if not feature_ids: break
        if features:
            output_indices = [i for i in output_indices if self.get_i(i).feature in features]
        ## prepare output
        output_indices = sorted(output_indices)
        if index: return output_indices
        else: return self.get_i(sorted(output_indices))
    
    def get_features_and_subfeatures(self, feature_ids, index = False, full = True):
        '''
        Gets features w/ feature_ids AND subfeatures of those features.
        If full = True, executes get_subfeatures_full for subfeature discovery, else get_subfeatures
        '''
        features = self.get_id(feature_ids, index = True, output_list = True)
        if full: subfeatures = self.get_subfeatures_full(feature_ids, index = True)
        else: subfeatures = self.get_subfeatures(feature_ids, index = True)
        ## prepare output
        final_features = sorted(features + subfeatures)
        if index: return final_features
        else: return self.get_i(final_features)

    def get_i_raw(self, indices, strip_newline = True):
        indices = indices if not isinstance(indices, int) else [indices]
        return self._indexed_file.get_line(*indices, strip_newline = strip_newline)
    
    def get_i(self, indices):
        '''
        Get GFF entry/entries by line index
        '''
        if self._data:
            if type(indices) is int: return self._data[indices]
            else: return [self._data[i] for i in indices]
        elif self.is_indexed():
            raw_entries = self.get_i_raw(indices, strip_newline = True)
            ## newline already stripped by get_i_raw. No need to strip them again.
            mk_annotation = self.make_annotation_from_str_gen(strip_newline = False)
            if type(indices) is int: return mk_annotation(raw_entries)
            else: return [mk_annotation(entry) for entry in raw_entries]
    
    def get_id(self, feature_ids, index = False, output_list = False):
        ## if even if output_list is only used when len(indices) == 1 or == 0 AND type(feature_ids) is str.
        ##  Always returns list otherwise.
        '''
        Get GFF entry by feature ID
        '''
        if type(feature_ids) is str:
            indices = [i for i, entry in enumerate(self) if entry.has_attr("ID", [feature_ids])]
        else:
            indices = [i for i, entry in enumerate(self) if entry.has_attr("ID", feature_ids)]
        if type(feature_ids) and not output_list and len(indices) <= 1:
            if not indices: return None
            elif index: return indices[0]
            else: return self.get_i(indices[0])
        elif index: return indices
        else: return list(self.get_i(indices))
    
    def write(self, fout, entries, **kwargs):
        '''
        Writes entries to file
        '''
        with open(fout, "w+") as f:
            f.write('\n'.join([entry.generate_str(**kwargs) for entry in entries]) + '\n')
    
    def write_i(self, fout, indices, **kwargs):
        '''
        Executes get_i, then writes output to file
        '''
        self.write(fout, self.get_i(indices), **kwargs)
        return
    
    def write_id(self, fout, feature_ids, **kwargs):
        '''
        Executes get_id, then writes output to file
        '''
        self.write(fout, self.get_id(feature_ids), **kwargs)
        
class Annotation:
    def __init__(self, entry, gff, **kwargs):
        self._gff = gff
        self.molecule = entry[0]
        self.source = entry[1]
        self.feature = entry[2]
        self.start = int(entry[3])
        self.end = int(entry[4])
        self.score = entry[5] if entry[5] == '.' else float(entry[5])
        self.strand = entry[6]
        self.phase = entry[7] if entry[7] == '.' else int(entry[7])
        self.attributes = Attributes(entry[8], gff, self, **kwargs)
        self._fields = {"molecule": self.molecule,
                        "source": self.source,
                        "feature": self.feature,
                        "start": self.start,
                        "end": self.end,
                        "score": self.score,
                        "strand": self.strand,
                        "phase": self.phase,
                        "attributes": self.attributes}
    def generate_attr(self, original = True, fields = None):
        if original: return self.attributes._raw
        else: return self.attributes.standardise_fields()
    def generate_str(self, fmt = "GFF"):
        if fmt.upper() in {"GFF", "GFF3"}:
            output = self.generate_gff()
        elif fmt.upper() in {"BED"}:
            output = self.generate_bed()
        return '\t'.join(map(str, output))
    def generate_gff(self):
        return list(map(str, [self.molecule, self.source, self.feature, self.start, self.end,
                              self.score, self.strand, self.phase, self.generate_attr()]))
    def generate_bed(self):
        return list(map(str, [self.molecule, self.start - 1, self.end, self.attributes.get("ID", fmt = str),
                              self.score, self.strand, self.source, self.feature, self.phase,
                              self.generate_attr()]))
    def get(self, *fields):
        return [self.f_dict[field] for field in fields]
    
    ## wrappers for Attribute methods
    def get_attr(self, a, **kwargs): return self.attributes.get(a, **kwargs)
    def is_attr(self, a, val, **kwargs): return self.attributes.is_attr(a, val, **kwargs)
    def has_attr(self, a, vals, **kwargs): return self.attributes.has_attr(a, vals, **kwargs)

class Attributes:
    def __init__(self, val, gff, entry, field_sep_inter = ';', field_sep_intra = ',',
                 **for_dummy_gff):
        self._entry = entry
        self._gff = gff if gff is not None else GFF(**for_dummy_gff)
        self._raw = val
        self._sep_inter = field_sep_inter
        self._sep_intra = field_sep_intra
        self._data = self._parse()
    def __repr__(self):
        return self._raw
    def __str__(self):
        return self._raw
    def _parse(self):
        import re
        ## if multiple separate entries for same field (e.g. 'Parent=abc.1;Parent=abc.2'), parse properly
        attributes_l = [x.split('=') for x in re.findall(f"[^{self._sep_intra}{self._sep_inter}=]+=.+?(?=[^{self._sep_intra}{self._sep_inter}=]+=.+?|$)", self._raw)]
        attributes = {}
        for attribute in attributes_l:
            attributes[attribute[0]] = attributes.get(attribute[0], []) + \
                                       re.search(f'^.+?(?=[{self._sep_intra}{self._sep_inter}]?$)',
                                                 attribute[1]).group(0).split(self._sep_intra)
        return attributes
    def get(self, a, fmt = list):
        feature = self._entry.feature
        attr_fields = self._gff._attr_fields
        if feature in attr_fields and a in attr_fields[feature]:
            mapped_field = attr_fields[feature][a]
        elif a in attr_fields["all"]:
            mapped_field = attr_fields["all"][a]
        else:
            mapped_field = a
        iterable = self._data.get(mapped_field, [])
        ## format return
        if fmt is str and iterable: return self._sep_intra.join(iterable)
        elif fmt is str and not iterable: return '.'
        else: return fmt(iterable)
    def is_attr(self, a, val):
        '''
        Checks if 'val' is a value of attribute 'a'
        '''
        return val in self.get(a)        
    def has_attr(self, a, vals):
        '''
        Checks if at least 1 value in 'vals' is also a value of attribute 'a'
        '''
        return (set(self.get(a)) & set(vals)) ## checks intersection != empty set
    def standardise_fields(self):
        def get_mod(field):
            feature_mod = get_recursively(self._gff._attr_fields_inv, field, self.feature, field)
            if feature_mod != field: return feature_mod
            else: return get_recursively(self._gff._attr_fields_inv, field, "all", field)
        return ';'.join([f"{get_mod(k)}={'.'.join(v)}" for k, v in self._data.items()])

# class GFF:
    
#     def __init__(self, fname = None, data = [], attr_mod = {}, fmt = "GFF3", **kwargs):
#         self._fmt = fmt
#         self._fname = fname
#         self._data = data
#         self._attr_mod = attr_mod
#         self._attr_fields = {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
#                                      "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
#                                      "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
#                                      "Is_circular": "Is_circular"}}
#         self._attr_fields_inv = {}
#         self._kwargs = kwargs
#         self.update_attr_fields()
#         self.parse()
    
#     def parse(self):
#         if self._fname is None: return
#         self._data = []
#         raw_data = [x.split('\t') for x in splitlines(self._fname) if x[:1] != '#']
#         if self._fmt.upper() == "BED":
#             def add_entry(entry):
#                 gff_fmt = [entry[0], entry[6], entry[7], str(int(entry[1]) + 1), entry[2],
#                            entry[4], entry[5], entry[8], entry[9]]
#                 self.add_entry(Annotation(gff_fmt, self, **self._kwargs))
#         else: ## GFF3
#             def add_entry(entry):
#                 self.add_entry(Annotation(entry, self, **self._kwargs))
#         for entry in [x for x in raw_data]:
#             add_entry(entry)
    
#     # def parse(self):
#     #     if self._fname is None: return
#     #     self._data = []
#     #     raw_data = [x.split('\t') for x in splitlines(self._fname) if x[:1] != '#']
#     #     for entry in [x for x in raw_data]:
#     #         self.add_entry(Annotation(entry, self, **self._kwargs))
    
#     def read_file(self, fname):
#         self._fname = fname
#         self.parse()
    
#     def add_entry(self, gff_entry):
#         self._data.append(gff_entry)
    
#     def update_attr_fields(self, attr_mod = None):
#         '''
#         Update attribute fields w/ user-provided attribute field modification dictionary.
#         (Otherwise, use self._attr_fields)
#         '''
#         if attr_mod is None: attr_mod = self._attr_mod
#         for feature, mods in attr_mod.items():
#             if feature in self._attr_fields:
#                 for canon, atypical in mods.items():
#                     self._attr_fields[feature][canon] = atypical
#             else:
#                 self._attr_fields[feature] = mods
#         self._attr_fields_inv = self.invert_attr_fields()
#         return
    
#     def invert_attr_fields(self):
#         output = {feature: {v: k for k, v in attr_map.items()}
#                  for feature, attr_map in self._attr_fields.items()}
#         return output
    
#     def get_subfeatures(self, feature_ids, *features, index = False):
#         '''
#         Get subfeatures of features w/ user-provided feature_ids.
#         '''
#         features = set(features)
#         feature_ids = set(feature_ids)
#         indices = [i for i, entry in enumerate(self._data) if
#                    ((not features or entry.feature in features) and
#                     entry.has_attr("Parent", feature_ids))]
#         if index: return indices
#         else: return [self._data[i] for i in indices]
    
#     def get_subfeatures_full(self, feature_ids, *features, index = False):
#         '''
#         Get all features that are subfeatures of user-provided feature_ids
#         AND subfeatures of those subfeatures, until there are no sub-sub...sub-features left.
#         '''
#         features = set(features)
#         output_indices = []
#         curr_indices = []
#         iteration_n = 0
#         while True:
#             iteration_n += 1
#             print(f"Executing iteration {iteration_n}")
#             curr_indices = self.get_subfeatures(feature_ids, index = True)
#             feature_ids = set(itertools.chain(*[self.get_i(i).get_attr("ID") for i in curr_indices]))
#             output_indices.extend(curr_indices)
#             if not feature_ids: break
#         if features:
#             output_indices = [i for i in output_indices if self.get_i(i).feature in features]
#         ## prepare output
#         output_indices = sorted(output_indices)
#         if index: return output_indices
#         else: return self.get_i(sorted(output_indices))
    
#     def get_features_and_subfeatures(self, feature_ids, index = False, full = True):
#         '''
#         Gets features w/ feature_ids AND subfeatures of those features.
#         If full = True, executes get_subfeatures_full for subfeature discovery, else get_subfeatures
#         '''
#         features = self.get_id(feature_ids, index = True)
#         if full: subfeatures = self.get_subfeatures_full(feature_ids, index = True)
#         else: subfeatures = self.get_subfeatures(feature_ids, index = True)
#         ## prepare output
#         final_features = sorted(features + subfeatures)
#         if index: return final_features
#         else: return self.get_i(final_features)
    
#     def get_i(self, indices):
#         '''
#         Get GFF entry by index
#         '''
#         if type(indices) is int: return self._data[indices]
#         else: return [self._data[i] for i in indices]
    
#     def get_id(self, feature_ids, index = False):
#         '''
#         Get GFF entry by feature ID
#         '''
#         indices = [i for i, entry in enumerate(self._data) if entry.has_attr("ID", feature_ids)]
#         if not indices and type(feature_ids) is str: return None
#         elif type(feature_ids) is str: return self.get_i(indices[0])
#         elif index: return indices
#         else: return [self.get_i(i) for i in indices]
    
#     def write(self, fout, entries, **kwargs):
#         '''
#         Writes entries to file
#         '''
#         open(fout, "w+").write('\n'.join([entry.generate_str(**kwargs) for entry in entries]) + '\n')
    
#     def write_i(self, fout, indices, **kwargs):
#         '''
#         Executes get_i, then writes output to file
#         '''
#         self.write(fout, self.get_i(indices), **kwargs)
#         return
    
#     def write_id(self, fout, feature_ids, **kwargs):
#         '''
#         Executes get_id, then writes output to file
#         '''
#         self.write(fout, self.get_id(feature_ids), **kwargs)
        


# class Annotation:
#     def __init__(self, entry, gff, **kwargs):
#         self._gff = gff
#         self.molecule = entry[0]
#         self.source = entry[1]
#         self.feature = entry[2]
#         self.start = int(entry[3])
#         self.end = int(entry[4])
#         self.score = entry[5] if entry[5] == '.' else float(entry[5])
#         self.strand = entry[6]
#         self.phase = entry[7] if entry[7] == '.' else int(entry[7])
#         self.attributes = Attributes(entry[8], gff, self, **kwargs)
#         self._fields = {"molecule": self.molecule,
#                         "source": self.source,
#                         "feature": self.feature,
#                         "start": self.start,
#                         "end": self.end,
#                         "score": self.score,
#                         "strand": self.strand,
#                         "phase": self.phase,
#                         "attributes": self.attributes}
#     def generate_attr(self, original = True, fields = None):
#         if original: return self.attributes._raw
#         else: return self.attributes.standardise_fields()
#     def generate_str(self, fmt = "GFF"):
#         if fmt.upper() in {"GFF", "GFF3"}:
#             output = self.generate_gff()
#         elif fmt.upper() in {"BED"}:
#             output = self.generate_bed()
#         return '\t'.join(map(str, output))
#     def generate_gff(self):
#         return list(map(str, [self.molecule, self.source, self.feature, self.start, self.end,
#                               self.score, self.strand, self.phase, self.generate_attr()]))
#     def generate_bed(self):
#         return list(map(str, [self.molecule, self.start - 1, self.end, self.attributes.get("ID", fmt = str),
#                               self.score, self.strand, self.source, self.feature, self.phase,
#                               self.generate_attr()]))
#     def get(self, *fields):
#         return [self.f_dict[field] for field in fields]
    
#     ## wrappers for Attribute methods
#     def get_attr(self, a): return self.attributes.get(a)
#     def is_attr(self, a, val): return self.attributes.is_attr(a, val)
#     def has_attr(self, a, vals): return self.attributes.has_attr(a, vals)

# class Attributes:
#     def __init__(self, val, gff, entry, field_sep_inter = ';', field_sep_intra = ','):
#         self._entry = entry
#         self._gff = gff
#         self._raw = val
#         self._sep_inter = field_sep_inter
#         self._sep_intra = field_sep_intra
#         self._data = self._parse()
#     def __repr__(self):
#         return self._raw
#     def __str__(self):
#         return self._raw
#     def _parse(self):
#         ## if multiple separate entries for same field (e.g. 'Parent=abc.1;Parent=abc.2'), parse properly
#         attributes_l = [x.split('=') for x in re.findall(f"[^{self._sep_intra}{self._sep_inter}=]+=.+?(?=[^{self._sep_intra}{self._sep_inter}=]+=.+?|$)", self._raw)]
#         attributes = {}
#         for attribute in attributes_l:
#             attributes[attribute[0]] = attributes.get(attribute[0], []) + \
#                                        re.search(f'^.+?(?=[{self._sep_intra}{self._sep_inter}]?$)',
#                                                  attribute[1]).group(0).split(self._sep_intra)
#         return attributes
#     def get(self, a, fmt = list):
#         feature = self._entry.feature
#         attr_fields = self._gff._attr_fields
#         if feature in attr_fields and a in attr_fields[feature]:
#             mapped_field = attr_fields[feature][a]
#         elif a in attr_fields["all"]:
#             mapped_field = attr_fields["all"][a]
#         else:
#             mapped_field = a
#         iterable = self._data.get(mapped_field, [])
#         ## format return
#         if fmt is str and iterable: return self._sep_intra.join(iterable)
#         elif fmt is str and not iterable: return '.'
#         else: return fmt(iterable)
#     def is_attr(self, a, val):
#         '''
#         Checks if 'val' is a value of attribute 'a'
#         '''
#         return val in self.get(a)
        
#     def has_attr(self, a, vals):
#         '''
#         Checks if at least 1 value in 'vals' is also a value of attribute 'a'
#         '''
#         return (set(self.get(a)) & set(vals)) ## checks intersection != empty set
#     def standardise_fields(self):
#         def get_mod(field):
#             feature_mod = get_recursively(self._gff._attr_fields_inv, field, self.feature, field)
#             if feature_mod != field: return feature_mod
#             else: return get_recursively(self._gff._attr_fields_inv, field, "all", field)
#         return ';'.join([f"{get_mod(k)}={'.'.join(v)}" for k, v in self._data.items()])
