import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from reference import AnnotatedFasta, ranges_union
from annotation import GFF, Annotation
from fasta_manip import dict_to_fasta
from pathlib import Path

# dir_results = "/mnt/chaelab/rachelle/nlr_survey/results"
# dir_out = dir_results + "/gff_curated"

def get_flanked_seq(ref, fids, flank = 200, feature_type = "gene",
                    gff_subset = None,
                    fasta_out = None, gff_out = None,
                    fasta_dir = None, gff_dir = None, gb_dir = None,
                    prefix = "flanked", adjust_dir = False, adjust_dir_fasta = False,
                    gff_supp = [], id_as_label = False, linked_features = [], unlinked_features = []):
    """
    Subsets GFF3 annotation by feature IDs, includes subfeatures.
    Extends feature range by 'flank' bp and subtracts the lower value of the extended range
    from all coordinates in feature itself + its subfeatures.
    Extracts extended sequences from FASTA file, named by feature ID.
    Writes gff3 subsetted from original to file (coordinates per original GFF3 file) if gff_subset != None.
    Writes gff3 w/ sequences (coordinates per extracted sequence) to file if gff_out != None.
    Writes sequences to file if fasta_out != None.
    Writes gff3 w/ sequence BY FEATURE ID if gff_dir != None (into {gff_dir}/{prefix}_{fid}.gff3).
    (Includes annotations of subfeatures)
    Writes FASTA format BY FEATURE ID if fasta_dir != None (into {fasta_dir}/{prefix}_{fid}.gff3).
    Writes GenBank format BY FEATURE ID if gb_dir != None (into {gb_dir}/{prefix}_{fid}.gff3).
    
    Arguments:
        ref (AnnotatedFasta): AnnotatedFasta object
        fids (list/tuple): reusable iterable of feature IDs
        flank (int): number of bp to extend range
        feature_type (str): GFF3 feature type (default="gene")
        gff_subset (str): path to output subsetted GFF3 file (coordinates per original GFF3 file) (optional)
        fasta_out (str): path to output FASTA file (write all subsetted sequences written to this single file) (optional)
        gff_out (str): path to output GFF3 file (coordinates per extracted sequence) (optional)
        fasta_dir (str): path to output FASTA directory (optional)
        gb_dir (str): path to output genbank directory (optional)
        prefix (str): prefix for files generated if gff_dir, fasta_dir, and/or gb_dir != None
        adjust_dir (bool): set to True to reverse complement sequence for output of fasta_out, fasta_dir, and gb_dir files AND adjust direction and strand in gff_out
        adjust_dir_fasta (bool): set to True to write plus strand of genes
        gff_supp (list): list of GFF class objects with additional features to be included if within range
            (for adding annotations like modified bases)
    
    Returns:
        dict: dictionary of extended sequences (format: {<feature ID>: <extended seq>})
    """
    ## subset original GFF3 file
    ref.subset_annotation(*fids)
    if gff_subset is not None:
        ref.annotation.write(fout = gff_subset)
    ## subset entries in gff_supp and add to reference annotation
    supp_gff = GFF()
    if gff_supp:
        ## get ranges to subset within
        ranges = {}
        for entry in ref.annotation:
            ranges[entry.molecule] = ranges.get(entry.molecule, []) + [(entry.start-flank, entry.end+flank)]
        ranges = {chrom: ranges_union([r]) for chrom, r in ranges.items()}
        ## subset
        for gff in gff_supp:
            supp_gff._data.extend(gff.get_in_range(ranges, output_list = True, index = False))
    ## extract sequences +
    ## generate GFF3 annotations with coordinates for extracted sequence
    new_gff = GFF()
    sub_gffs = {}
    seq_output = {}
    seq_output_adj = {}
    for fid in fids:
        entry = ref.annotated_feature(fid, feature_type = feature_type)
        molecule, start, end = entry.molecule, entry.start, entry.end
        true_start = max(0, start-flank)
        seq = ref[molecule][true_start:end+flank]
        seq_output[fid] = seq
        if adjust_dir or ((fasta_out or fasta_dir) and adjust_dir_fasta):
            seq_output_adj[fid] = seq.reverse_complement() if entry.strand == '-' else seq
        ## subset main GFF entries
        entries = ref.annotation.subset(fid, subfeatures = True)
        sub_gff = GFF()
        for entry in entries:
            ## create new Annotation object
            new_entry = Annotation(entry.generate_str(fmt = "GFF3").split('\t'), new_gff)
            new_entry.start = entry.start - true_start
            new_entry.end = entry.end - true_start
            if adjust_dir and entry.strand == '-':
                new_start = new_entry.start
                new_end = new_entry.end
                new_entry.start = len(seq) - new_end
                new_entry.end = len(seq) - new_start
                new_entry.strand = '+'
            new_entry.seqid = fid
            sub_gff._data.append(new_entry)
            new_gff._data.append(new_entry)
        ## subset supp_gff entries
        supp_entries = supp_gff.subset(ranges = {molecule: [(true_start, end+flank)]})
        for entry in supp_entries:
            ## create new Annotation object
            new_entry = Annotation(entry.generate_str(fmt = "GFF3").split('\t'), new_gff)
            new_entry.start = entry.start - true_start
            new_entry.end = entry.end - true_start
            if adjust_dir and entry.strand == '-':
                new_start = new_entry.start
                new_end = new_entry.end
                new_entry.start = len(seq) - new_end
                new_entry.end = len(seq) - new_start
                new_entry.strand = '+'
            new_entry.seqid = fid
            new_entry.attributes._data["Parent"] = [fid] ## this is required to link the annotation to the feature
            new_entry.attributes._raw = new_entry.generate_str(standardise = True)
            sub_gff._data.append(new_entry)
            new_gff._data.append(new_entry)
            ## add original entry to ref.annotation
            ref.annotation._data.append(entry)
        ## add to sub_gffs
        sub_gffs[fid] = sub_gff
    ## write again, this time updated with gff_supp entries
    if gff_subset is not None:
        ref.annotation.write(fout = gff_subset)
    if gff_out is not None:
        new_gff.write(fout = gff_out, seqs = (seq_output_adj if adjust_dir else seq_output))
    if gff_dir is not None:
        Path(gff_dir).mkdir(parents = True, exist_ok = True)
        for fid, sub_gff in sub_gffs.items():
            sub_gff.write(fout = f"{gff_dir}/{prefix}_{fid}.gff3",
                          seqs = {fid: (seq_output_adj if adjust_dir else seq_output)[fid]})
    seq_output_for_fasta = seq_output_adj if (adjust_dir or adjust_dir_fasta) else seq_output
    if fasta_out is not None:
        dict_to_fasta(seq_output_for_fasta, fasta_out)
    if fasta_dir is not None:
        Path(fasta_dir).mkdir(parents = True, exist_ok = True)
        for seqid, seq in seq_output_for_fasta.items():
            dict_to_fasta({seqid: seq}, f"{fasta_dir}/{prefix}_{seqid}.fasta")
    if gb_dir is not None:
        Path(gb_dir).mkdir(parents = True, exist_ok = True)
        for fid, sub_gff in sub_gffs.items():
            print(fid)
            with open(f"{gb_dir}/{prefix}_{fid}.gb", 'w+') as f:
                f.write(make_genbank_for_gene(fid, (seq_output_adj if adjust_dir else seq_output)[fid], sub_gff,
                                              id_as_label = id_as_label,
                                              linked_features = linked_features,
                                              unlinked_features = unlinked_features))
    return seq_output

# seqs_flanked = get_flanked_seq(bra3, bra_rnl_id, flank = 200, feature_type = "gene",
#                                gff_out = f"{dir_out}/Brapa3_ADR1_NRG1.gff3",
#                                gff_dir = f"{dir_out}/separate",
#                                fasta_out = f"{dir_out}/Brapa3_ADR1_NRG1.fasta",
#                                fasta_dir = f"{dir_out}/separate")

## takes list of Annotation objects and makes a GenBank feature
## mostly to be used for generating string output tbh
class GenBankFeature:
    def __init__(self, gene, type, annotations, parent = None, seq = '', **kwargs):
        self.type = type
        self.gene = gene ## str of ID
        self.parent = parent ## str of ID
        self.annotations = annotations
        self._seq = seq
        self.data = kwargs
        self.plus = self.annotations[0].strand == '+'
    def ranges(self):
        return sorted([(e.start, e.end) for e in self.annotations])
    def seq(self):
        ranges = self.ranges()
        output = self._seq[ranges[0][0]:ranges[0][1]]
        for r in ranges[1:]:
            output += self._seq[r[0]:r[1]]
        return output
    def translate(self):
        nuc = self.seq()
        return nuc.translate() if self.plus else nuc.reverse_complement().translate()
    def format_range(self):
        ranges = sorted(self.ranges(), reverse = (not self.plus))
        sranges = ([f"{r[0]+1}..{r[1]}" for r in ranges] if self.plus
                   else [f"complement({r[0]+1}..{r[1]})" for r in ranges])
        return (sranges[0] if len(sranges) == 1
                else f"join({','.join(sranges)})")
    def format_entry(self, id_as_label = False):
        lines = []
        lines.append( (' '*5) + self.type + (' '*(16-len(self.type))) + self.format_range())
        def make_line(key, value):
            lines.append((' '*21) + '/' + key + '=' +
                         (f'"{value}"' if isinstance(value, str) else str(value)))
        make_line("gene", self.gene)
        if self.type == "CDS":
            parents = self.annotations[0].get_attr("Parent")
            make_line("codon_start",
                      (min(self.annotations, key = lambda a:a.start) if self.plus
                       else max(self.annotations, key = lambda a:a.end)).phase)
            make_line("product", parents[0])
            if id_as_label: make_line("label", parents[0])
            if self._seq: make_line("translation", self.translate())
        else:
            ids = self.annotations[0].get_attr("ID", fmt = list)
            if self.type == "mRNA":
                if ids: make_line("product", ids[0])
            if id_as_label and ids: make_line("label", ids[0])
        return '\n'.join(lines)
    def format_seq(self):
        output = ''
        for i in range(0, len(self._seq), 60):
            line = ' '*(9-len(str(i))) + str(i)
            for j in range(i, i+60, 10):
                line += ' ' + str(self._seq[j:j+10])
            output += line + '\n'
        return output[:-1]

## only parses the following feature types: gene, mRNA, CDS
## also only works for EXACTLY ONE gene (1 gid). Any more or less, and the output may be unexpected
## linked_features and unlinked_features are for any other feature that should be included (e.g. UTR)
##  linked_features is for features that should be joined (e.g. something like CDS) (in beta, no way to separate different features of the same type for now)
##  unlinked_features is for features that are disjointed (e.g. modified bases)
#### ignore ##  format: {<parent feature type>: [<feature type 1>, <feature type 2>]}
def make_genbank_for_gene(gid, seq, gff, unlinked_features = [], linked_features = [], id_as_label = False):
    output = []
    gff = gff.subset(feature_ids = [gid], subfeatures = True)
    ## header
    import datetime
    output.append("LOCUS" + ' '*7 + gid + (' '*(28-len(gid)-len(str(len(seq))))) +
                  str(len(seq)) + " bp" + "     DNA" + "    linear" + "   UNC" +
                  datetime.datetime.now().strftime("%d-%b-%Y").upper())
    ## format features
    output.append("FEATURES")
    gene = [a for a in gff if a.type == "gene"]
    genbank_gene = GenBankFeature(gid, "gene", gene, parent = None, seq = seq)
    output.append(genbank_gene.format_entry(id_as_label = id_as_label))
    mrnas = [a for a in gff if a.type == "mRNA"]
    for mrna in mrnas:
        output.append(GenBankFeature(gid, "mRNA", [mrna], parent = gid, seq = seq).format_entry(id_as_label = id_as_label))
        cds = gff.get_subfeatures(*mrna.get_attr("ID"), feature_types = ["CDS"])
        if cds:
            output.append(GenBankFeature(gid, "CDS", cds,
                                         parent = mrna.get_attr("ID")[0], seq = seq).format_entry(id_as_label = id_as_label))
    for feature_type in unlinked_features:
        to_add = gff.get_features(feature_type, index = False)
        if len(to_add) > 0:
            for entry in to_add:
                output.append(GenBankFeature(gid, feature_type, [entry], parent = gid, seq = seq).format_entry(id_as_label = id_as_label))
    for feature_type in linked_features:
        to_add = gff.get_features(feature_type, index = False)
        if len(to_add) > 0:
            output.append(GenBankFeature(gid, feature_type, to_add, parent = gid, seq = seq).format_entry(id_as_label = id_as_label))
    ## format sequence
    output.append("ORIGIN")
    output.append(genbank_gene.format_seq())
    output.append("\\\\")
    return '\n'.join(output)
        

# test = GFF(fname = "/mnt/chaelab/rachelle/nlr_survey/results/gff_curated/separate/flanked_BraA06g033140.3C.gff3")
# ftest = fasta_to_dict("/mnt/chaelab/rachelle/nlr_survey/results/gff_curated/separate/flanked_BraA06g033140.3C.fasta")
# # with open(dir_out + "/test/test.gb", 'w+') as f:
# #     f.write(make_genbank_for_gene("BraA06g033140.3C", ftest["BraA06g033140.3C"], test))

# dir_brapa = "/mnt/chaelab/rachelle/data/Brapa"

# bra3 = AnnotatedFasta(f"{dir_brapa}/v3.0/Brapa_sequence_v3.0.fasta",
#                       f"{dir_brapa}/v3.0/Brapa_genome_v3.0_genes.gff3",
#                       name = "Brapa3")
# ## gene ids extracted using the following script in directory /mnt/chaelab/rachelle/nlr_survey/results/id_curated
# ## cd /mnt/chaelab/rachelle/nlr_survey/results/id_curated; echo "[\"$(grep BraA Col0-Brapa3-Cai2021.apt_TCN.*.id.txt | tr ':' '\t' | cut -f2 | tr '\n' ',' | sed 's/,$//' | sed 's/,/","/g')\"]"
# bra_rnl_id = ["BraA01g004700.3C","BraA08g016670.3C","BraA10g031680.3C","BraA02g040960.3C","BraA06g033140.3C","BraA06g037470.3C","BraA07g017010.3C","BraA07g017020.3C","BraA09g009300.3C","BraA09g042720.3C"]
# # bra3.subset_annotation(*bra_rnl_id)

# seqs_flanked = get_flanked_seq(bra3, bra_rnl_id, flank = 200, feature_type = "gene",
#                                prefix = "Bra3_flank200bp",
#                                gff_subset = f"{dir_out}/Brapa3_ADR1_NRG1.gff3",
#                                # gff_out = f"{dir_out}/Brapa3_ADR1_NRG1.gff3",
#                                gff_dir = f"{dir_out}/RNL",
#                                # fasta_out = f"{dir_out}/Brapa3_ADR1_NRG1.fasta",
#                                # fasta_dir = f"{dir_out}/separate",
#                                gb_dir = f"{dir_results}/gb_curated/RNL")

# dir_tcnl = "/mnt/chaelab/rachelle/tgfam_finder/results/arath_potato_tomato_TCNLdomains"

# ## cd /mnt/chaelab/rachelle/nlr_survey/results/id_curated; echo "[\"$(grep R500 Col0-Brapa3-R500-Z1-Cai2021.apt_TCNL.*.id.txt | tr ':' '\t' | cut -f2 | tr '\n' ',' | sed 's/,$//' | sed 's/,/","/g')\"]"
# r500_rnl_id = ["Brapa_R500_NLR.A02.1","Brapa_R500_NLR.A08.7","Brapa_R500_NLR.A10.10","Brapa_R500_NLR.A02.40","Brapa_R500_NLR.A06.14","Brapa_R500_NLR.A06.16","Brapa_R500_NLR.A07.2","Brapa_R500_NLR.A07.3","Brapa_R500_NLR.A09.2"]
# z1_rnl_id = ["Brapa_Z1_NLR.A01.3","Brapa_Z1_NLR.A08.8","Brapa_Z1_NLR.A10.11","Brapa_Z1_NLR.A02.35","Brapa_Z1_NLR.A06.18","Brapa_Z1_NLR.A06.23","Brapa_Z1_NLR.A06.24","Brapa_Z1_NLR.A07.2","Brapa_Z1_NLR.A07.3","Brapa_Z1_NLR.A09.3"]
# r500 = AnnotatedFasta(f"{dir_brapa}/R500/Brapa_R500_v1.2.fasta",
#                       f"{dir_tcnl}/Brapa_R500_TCNL/Final/Brapa_R500_NLR.TGFam.hierarchy.gff3",
#                       name = "R500")
# z1 = AnnotatedFasta(f"{dir_brapa}/Z1/GCA_900412535.2_Brassica_rapa_Z1_genomic_renamed.fna",
#                     f"{dir_tcnl}/Brapa_Z1_TCNL/Final/Brapa_Z1_NLR.TGFam.hierarchy.gff3",
#                     name = "R500")

# seqs_flanked = get_flanked_seq(r500, r500_rnl_id, flank = 200, feature_type = "gene",
#                                prefix = "R500_tgfam-TCNL_flank200bp",
#                                gff_subset = f"{dir_out}/R500_tgfam-TCNL_ADR1_NRG1.gff3",
#                                # gff_out = f"{dir_out}/R500_tgfam-TCNL_ADR1_NRG1.gff3",
#                                gff_dir = f"{dir_out}/RNL",
#                                gb_dir = f"{dir_results}/gb_curated/RNL")
# seqs_flanked = get_flanked_seq(z1, z1_rnl_id, flank = 200, feature_type = "gene",
#                                prefix = "Z1_tgfam-TCNL_flank200bp",
#                                gff_subset = f"{dir_out}/Z1_tgfam-TCNL_ADR1_NRG1.gff3",
#                                # gff_out = f"{dir_out}/Z1_tgfam-TCNL_ADR1_NRG1.gff3",
#                                gff_dir = f"{dir_out}/RNL",
#                                gb_dir = f"{dir_results}/gb_curated/RNL")


# dir_tcrnl = "/mnt/chaelab/rachelle/tgfam_finder/results/arath_potato_tomato_TCRNLdomains"

# ## cd /mnt/chaelab/rachelle/nlr_survey/results/id_curated; echo "[\"$(grep R500 Col0-Brapa3-R500-Z1-Cai2021.apt_TCRNL.*.id.txt | tr ':' '\t' | cut -f2 | tr '\n' ',' | sed 's/,$//' | sed 's/,/","/g')\"]"
# r500_rnl_id = ["Brapa_R500_NLR.A02.1","Brapa_R500_NLR.A08.7","Brapa_R500_NLR.A10.10","Brapa_R500_NLR.A02.40","Brapa_R500_NLR.A06.14","Brapa_R500_NLR.A06.16","Brapa_R500_NLR.A07.2","Brapa_R500_NLR.A07.3","Brapa_R500_NLR.A09.2","Brapa_R500_NLR.A09.54","Brapa_R500_NLR.A09.58"]
# z1_rnl_id = ["Brapa_Z1_NLR.A01.3","Brapa_Z1_NLR.A08.8","Brapa_Z1_NLR.A10.11","Brapa_Z1_NLR.A02.20","Brapa_Z1_NLR.A02.36","Brapa_Z1_NLR.A06.18","Brapa_Z1_NLR.A06.23","Brapa_Z1_NLR.A06.24","Brapa_Z1_NLR.A07.2","Brapa_Z1_NLR.A07.3","Brapa_Z1_NLR.A09.3"]
# r500 = AnnotatedFasta(f"{dir_brapa}/R500/Brapa_R500_v1.2.fasta",
#                       f"{dir_tcrnl}/Brapa_R500_TCRNL/Final/Brapa_R500_NLR.TGFam.hierarchy.gff3",
#                       name = "R500")
# z1 = AnnotatedFasta(f"{dir_brapa}/Z1/GCA_900412535.2_Brassica_rapa_Z1_genomic_renamed.fna",
#                     f"{dir_tcrnl}/Brapa_Z1_TCRNL/Final/Brapa_Z1_NLR.TGFam.hierarchy.gff3",
#                     name = "R500")

# seqs_flanked = get_flanked_seq(r500, r500_rnl_id, flank = 200, feature_type = "gene",
#                                prefix = "R500_tgfam-TCRNL_flank200bp",
#                                gff_subset = f"{dir_out}/R500_tgfam-TCRNL_ADR1_NRG1.gff3",
#                                # gff_out = f"{dir_out}/R500_tgfam-TCRNL_ADR1_NRG1.gff3",
#                                gff_dir = f"{dir_out}/RNL",
#                                gb_dir = f"{dir_results}/gb_curated/RNL")
# seqs_flanked = get_flanked_seq(z1, z1_rnl_id, flank = 200, feature_type = "gene",
#                                prefix = "Z1_tgfam-TCRNL_flank200bp",
#                                gff_subset = f"{dir_out}/Z1_tgfam-TCRNL_ADR1_NRG1.gff3",
#                                # gff_out = f"{dir_out}/Z1_tgfam-TCRNL_ADR1_NRG1.gff3",
#                                gff_dir = f"{dir_out}/RNL",
#                                gb_dir = f"{dir_results}/gb_curated/RNL")


# def get_flanked_seq_merge(ref, fids, flank = 200, feature_type = "gene",
#                           gff_subset = None,
#                           fasta_out = None, gff_out = None,
#                           fasta_dir = None, gff_dir = None, gb_dir = None,
#                           prefix = "flanked", adjust_dir_fasta = False):
#     """
#     Subsets GFF3 annotation by feature IDs, includes subfeatures.
#     Extends feature range by 'flank' bp and merges all overlapping or bookended sequences.
#     Subtracts the lower value of the extended range
#     from all coordinates in feature itself + its subfeatures.
#     Extracts extended sequences from FASTA file, named by feature ID.
#     Writes gff3 subsetted from original to file (coordinates per original GFF3 file) if gff_subset != None.
#     Writes gff3 w/ sequences (coordinates per extracted sequence) to file if gff_out != None.
#     Writes sequences to file if fasta_out != None.
#     Writes gff3 w/ sequence BY FEATURE ID if gff_dir != None (into {gff_dir}/{prefix}_{fid}.gff3).
#     (Includes annotations of subfeatures)
#     Writes FASTA format BY FEATURE ID if fasta_dir != None (into {fasta_dir}/{prefix}_{fid}.gff3).
#     Writes GenBank format BY FEATURE ID if gb_dir != None (into {gb_dir}/{prefix}_{fid}.gff3).
    
#     Arguments:
#         ref (AnnotatedFasta): AnnotatedFasta object
#         fids (list/tuple): reusable iterable of feature IDs
#         flank (int): number of bp to extend range
#         feature_type (str): GFF3 feature type (default="gene")
#         gff_subset (str): path to output subsetted GFF3 file (coordinates per original GFF3 file) (optional)
#         gff_out (str): path to output GFF3 file (coordinates per extracted sequence) (optional)
#         fasta_dir (str): path to output FASTA directory (optional)
#         gb_dir (str): path to output genbank directory (optional)
#         prefix (str): prefix for files generated if gff_dir, fasta_dir, and/or gb_dir != None
#         adjust_dir_fasta (bool): set to True to write plus strand of genes
    
#     Returns:
#         dict: dictionary of extended sequences (format: {<feature ID>: <extended seq>})
#     """
#     ## subset original GFF3 file
#     ref.subset_annotation(*fids)
#     if gff_subset is not None:
#         ref.annotation.write(fout = gff_subset)
#     ## extract sequences +
#     ## generate GFF3 annotations with coordinates for extracted sequence
#     new_gff = GFF()
#     sub_gffs = {}
#     seq_output = {}
#     seq_output_adj = {}
#     for fid in fids:
#         entry = ref.annotated_feature(fid, feature_type = feature_type)
#         mol, start, end = entry.molecule, entry.start, entry.end
#         true_start = max(0, start-flank)
#         seq = ref[mol][true_start:end+flank]
#         seq_output[fid] = seq
#         if (fasta_out or fasta_dir) and adjust_dir_fasta:
#             seq_output_adj[fid] = seq.reverse_complement() if entry.strand == '-' else seq
#         entries = ref.annotation.subset(fid)
#         sub_gff = GFF()
#         for entry in entries:
#             ## create new Annotation object
#             new_entry = Annotation(entry.generate_str(fmt = "GFF3").split('\t'), new_gff)
#             new_entry.start = entry.start - true_start
#             new_entry.end = entry.end - true_start
#             new_entry.seqid = fid
#             sub_gff._data.append(new_entry)
#             new_gff._data.append(new_entry)
#         sub_gffs[fid] = sub_gff
#     if gff_out is not None:
#         new_gff.write(fout = gff_out, seqs = seq_output)
#     if gff_dir is not None:
#         Path(gff_dir).mkdir(parents = True, exist_ok = True)
#         for fid, sub_gff in sub_gffs.items():
#             sub_gff.write(fout = f"{gff_dir}/{prefix}_{fid}.gff3",
#                           seqs = {fid: seq_output[fid]})
#     seq_output_for_fasta = seq_output if not adjust_dir_fasta else seq_output_adj
#     if fasta_out is not None:
#         dict_to_fasta(seq_output_for_fasta, fasta_out)
#     if fasta_dir is not None:
#         Path(fasta_dir).mkdir(parents = True, exist_ok = True)
#         for seqid, seq in seq_output_for_fasta.items():
#             dict_to_fasta({seqid: seq}, f"{fasta_dir}/{prefix}_{seqid}.fasta")
#     if gb_dir is not None:
#         Path(gb_dir).mkdir(parents = True, exist_ok = True)
#         for fid, sub_gff in sub_gffs.items():
#             print(fid)
#             with open(f"{gb_dir}/{prefix}_{fid}.gb", 'w+') as f:
#                 f.write(make_genbank_for_gene(fid, seq_output[fid], sub_gff))
#     return seq_output
