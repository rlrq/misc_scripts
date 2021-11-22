import re
import csv
import sys
import itertools

sys.path.append("/mnt/chaelab/rachelle/src")
from basics import get_recursively
from data_manip import splitlines

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
    
    def __init__(self, fname = None, data = [], attr_mod = {}, fmt = "GFF3", **kwargs):
        self._fmt = fmt
        self._fname = fname
        self._data = data
        self._attr_mod = attr_mod
        self._attr_fields = {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                                     "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                                     "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                                     "Is_circular": "Is_circular"}}
        self._attr_fields_inv = {}
        self._kwargs = kwargs
        self.update_attr_fields()
        self.parse()
    
    def parse(self):
        if self._fname is None: return
        self._data = []
        raw_data = [x.split('\t') for x in splitlines(self._fname) if x[:1] != '#']
        if self._fmt.upper() == "BED":
            def add_entry(entry):
                gff_fmt = [entry[0], entry[6], entry[7], str(int(entry[1]) + 1), entry[2],
                           entry[4], entry[5], entry[8], entry[9]]
                self.add_entry(Annotation(gff_fmt, self, **self._kwargs))
        else: ## GFF3
            def add_entry(entry):
                self.add_entry(Annotation(entry, self, **self._kwargs))
        for entry in [x for x in raw_data]:
            add_entry(entry)
    
    # def parse(self):
    #     if self._fname is None: return
    #     self._data = []
    #     raw_data = [x.split('\t') for x in splitlines(self._fname) if x[:1] != '#']
    #     for entry in [x for x in raw_data]:
    #         self.add_entry(Annotation(entry, self, **self._kwargs))
    
    def read_file(self, fname):
        self._fname = fname
        self.parse()
    
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
        feature_ids = set(feature_ids)
        indices = [i for i, entry in enumerate(self._data) if
                   ((not features or entry.feature in features) and
                    entry.has_attr("Parent", feature_ids))]
        if index: return indices
        else: return [self._data[i] for i in indices]
    
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
        features = self.get_id(feature_ids, index = True)
        if full: subfeatures = self.get_subfeatures_full(feature_ids, index = True)
        else: subfeatures = self.get_subfeatures(feature_ids, index = True)
        ## prepare output
        final_features = sorted(features + subfeatures)
        if index: return final_features
        else: return self.get_i(final_features)
    
    def get_i(self, indices):
        '''
        Get GFF entry by index
        '''
        if type(indices) is int: return self._data[indices]
        else: return [self._data[i] for i in indices]
    
    def get_id(self, feature_ids, index = False):
        '''
        Get GFF entry by feature ID
        '''
        indices = [i for i, entry in enumerate(self._data) if entry.has_attr("ID", feature_ids)]
        if not indices and type(feature_ids) is str: return None
        elif type(feature_ids) is str: return self.get_i(indices[0])
        elif index: return indices
        else: return [self.get_i(i) for i in indices]
    
    def write(self, fout, entries, **kwargs):
        '''
        Writes entries to file
        '''
        open(fout, "w+").write('\n'.join([entry.generate_str(**kwargs) for entry in entries]) + '\n')
    
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
    def get_attr(self, a): return self.attributes.get(a)
    def is_attr(self, a, val): return self.attributes.is_attr(a, val)
    def has_attr(self, a, vals): return self.attributes.has_attr(a, vals)

class Attributes:
    def __init__(self, val, gff, entry, field_sep_inter = ';', field_sep_intra = ','):
        self._entry = entry
        self._gff = gff
        self._raw = val
        self._sep_inter = field_sep_inter
        self._sep_intra = field_sep_intra
        self._data = self._parse()
    def __repr__(self):
        return self._raw
    def __str__(self):
        return self._raw
    def _parse(self):
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
