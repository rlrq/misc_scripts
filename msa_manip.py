def gapped_dna_msa_from_dict(d):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio.Alphabet import IUPAC, Gapped
    return MultipleSeqAlignment([SeqRecord(Seq(v, Gapped(IUPAC.unambiguous_dna)), id = k, description = "")
                                 for k, v in d.items()])
