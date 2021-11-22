import re

# wkdir = "/mnt/chaelab/rachelle/variants"
# # snpEff_f = wkdir + "/data/762genomes_snpEffects_10.bed"
# snpEff_f = wkdir + "/results/762genomes_snpEffects.bed"
# fout = wkdir + "/results/762genomes_snpEffects_expanded.bed"

## separate each entry into multiple lines based on effects annotation (ANN, comma separated), and also split the subfields into tab separated values
## valid only for Arabidopsis naming conventions
## compatible with snpEff v4.3
def expand_snpEff(snpEff_f, fout, fout_lof = "", generate_lof_file = True):
    snpEff_raw = []
    snpEff_dat = ['\t'.join(["chrom","start","end","ref","alt","effect","impact","geneName","geneIDdef","featureType","featureID","transcriptBiotype","rankTotal","HGVS.c","HGVS.p","cDNAposLen","CDSposLen","proteinPosLen","distToFeat","Errors","geneID"])]
    snpEff_lof = ['\t'.join(["chrom","start","end","ref","alt","LOF.NMD"])]
    with open(snpEff_f,'r') as f:
        snpEff_raw = [x[:-1].split('\t') for x in f.readlines()]
    for i,line in enumerate(snpEff_raw):
        if i+1 % 100000 == 0:
            print("{} lines processed".format(i))
        anns = line[3].split(';')
        anns_ref = anns[0]
        anns_main = anns[1][4:].split(',')
        anns_lof = anns[2:]
        for ann in anns_main:
            gene = re.search("AT\dG\d+",ann).group(0)
            snpEff_dat.append('\t'.join(line[:-1] + [anns_ref] +  ann.split('|') + [gene]))
        if anns_lof:
            snpEff_lof.append('\t'.join(line[:-1] + [anns_ref, anns[0].split('|')[0]] + anns_lof))
    import os
    f = open(fout,"w+")
    f.write('\n'.join(snpEff_dat))
    f.close()
    if generate_lof_file:
        f = open(fout_lof,"w+")
        f.write('\n'.join(snpEff_lof))
        f.close()
    print("Successfully expanded the following file: {}\nOutput files: {}".format(snpEff_f, fout) + "" if not generate_lof_file else ", {}".format(fout_lof))
