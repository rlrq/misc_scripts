## make ID map for non-paired IDs
def make_unpaired_id_map(f_paired_ids, f_out):
    with open(f_paired_ids, "r") as fin:
        with open(f_out, "w+") as fout:
            for i, inpt in enumerate(fin):
                _ = fout.write(str(i) + inpt.split(':')[1])
    return

# dir_in = "/mnt/chaelab/rachelle/panNLRome/results/regroup/OrthoFinder/vdw/vdw-Araport11.NLR.pep/blast_1.5/Results_Feb12/WorkingDirectory"
# f_in = dir_in + "/SequenceIDs.txt"
# f_out = dir_in + "/SequenceIDs.unpaired.txt"
# make_unpaired_id_map(f_in, f_out)

