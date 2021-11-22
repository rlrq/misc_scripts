import re
import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from basics import *
from data_manip import *

gff_colnames = ("molecule", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff_get = make_custom_get(gff_colnames)

gffbed_colnames = ("molecule", "start", "end", "name", "score", "strand", "source", "type", "phase",
                   "attributes")
gffbed_get = make_custom_get(gffbed_colnames)

# gff_fname_rice = "/mnt/chaelab/rachelle/data/rice/rice_nlr495.gff"
# fout_rice = lambda x: f"/mnt/chaelab/rachelle/data/rice/rice_nlr495_cluster{x}kb.txt"
# assign_clusters(gff_fname_rice, fout_rice, 50)
# assign_clusters(gff_fname_rice, fout_rice, 300)

# gffbed_fname_ricecomb = "/mnt/chaelab/rachelle/data/rice/rice_nlr596_MSU_TGFamFinder.GFF3.bed"
# fout_ricecomb = lambda x: f"/mnt/chaelab/rachelle/data/rice/rice_nlr596_MSU_TGFamFinder_cluster{x}kb.txt"
# assign_clusters(gffbed_fname_ricecomb, fout_ricecomb, 50, inpt = "bed", version = 2)

# gffbed_fname_ricecombext = "/mnt/chaelab/rachelle/data/rice/rice_nlr596ext_MSU_TGFamFinder.GFF3.bed"
# fout_ricecombext = lambda x: f"/mnt/chaelab/rachelle/data/rice/rice_nlr596ext_MSU_TGFamFinder_cluster{x}kb.txt"
# assign_clusters(gffbed_fname_ricecombext, fout_ricecombext, 50, inpt = "bed", version = 2)

# gffbed_fname_ricecombext = "/mnt/chaelab/rachelle/data/rice/rice_nlr596ext_MSU_TGFamFinder.GFF3.bed"
# unann_pred = "/mnt/chaelab/rachelle/cluster_survey/results/rice/MSU_v7_TGFam_IRGSP-1.0/Cdd.v3.18/alignment/rice_nlr596ext_ref-PRJEB23459_RX_superfam_pruned_seqnames.txt"
# fout_ricecombextpred = lambda x: f"/mnt/chaelab/rachelle/data/rice/rice_nlr596ext_MSU_TGFamFinder_pred_cluster{x}kb.txt"
# assign_clusters(gffbed_fname_ricecombext, fout_ricecombext, 50, inpt = "bed", version = 2)

# gffbed_fname_braptirnaz = "/mnt/chaelab/rachelle/data/Brapa/brapa_tirnaz_v3.0.GFF3.bed"
# fout_braptirnaz = lambda x: f"/mnt/chaelab/rachelle/data/Brapa/brapa_tirnaz_v3.0_cluster{x}kb.txt"
# assign_clusters_from_gff(gffbed_fname_braptirnaz, inpt = "bed",
#                          fout = fout_braptirnaz, cluster_within_kb = 50, version = 2)

gffbed_fname_braptirnaz = "/mnt/chaelab/rachelle/data/Brapa/brapa_tirnazext_v3.0.GFF3.bed"
fout_braptirnazext = lambda x: f"/mnt/chaelab/rachelle/data/Brapa/brapa_tirnazext_v3.0_cluster{x}kb.txt"
assign_clusters_from_gff(gffbed_fname_braptirnaz, inpt = "bed",
                         fout = fout_braptirnazext, cluster_within_kb = 50, version = 2)


## assumes gene pattern is AT.G.... or Os.g.....
def assign_clusters_from_gff(gff_fname, inpt = "gff", **kwargs):
    
    ## parse data
    dat = [x.split('\t') for x in splitlines(gff_fname)]
    dat_get = gff_get if inpt == "gff" else gffbed_get
    output = [list(dat_get(x, "molecule", "start", "end", "strand")) + \
              # [re.search("(?<=ID=).+?(?=;)", dat_get(x, "attributes")).group(0)] for x in dat
              [re.search("ID=([^;]+)", dat_get(x, "attributes")).group(1)] for x in dat
              if dat_get(x, "type") == "gene"]
    output.sort()
    # print(head(output))
    
    assign_clusters(output, **kwargs)
    return

def assign_clusters(dat, fout, cluster_within_kb, version = 1):
    
    cluster_within = cluster_within_kb * 1000
    colnames = ("chrom", "start", "end", "strand", "gene", "cluster")
    get = make_custom_get(colnames)
    dat = [get(x, "chrom", "start", "end", "strand", "gene") for x in dat]
    dat.sort()
    
    ## function to make and assign cluster name
    def assign_cluster(start, end, clust_num, # prefix = '', chrprefix = '[cC]hr',
                       version = 1): ## end is non-inclusive
        if version == 1:
            ## cluster name pattern: Os01c05
            clust_name = re.match(".+?(?=[gG])", get(dat[i], "gene")).group(0) + 'c' + str(clust_num).zfill(2)
        else:
            clust_name = f"{get(dat[i], 'chrom')}c{str(clust_num).zfill(2)}" ## just prepend 'c'
        for j in range(start, end):
            dat[j].append(clust_name)
    
    ## instantiate some variables
    last_entry = dat[0]
    clust_start = 0
    clust_num = 1
    n_entries = len(dat)
    
    ## start grouping genes into clusters
    for i, entry in enumerate(dat):
        is_last = i == (n_entries - 1)
        same_chrom = get(entry, "chrom") == get(last_entry, "chrom")
        within_range = get(entry, "start") - get(last_entry, "end") < cluster_within
        if not (same_chrom and within_range): ## different cluster
            if i - clust_start == 1:
                dat[clust_start].append("singleton")
                ## reset cluster number
                if not same_chrom:
                    clust_num = 1
            else:
                assign_cluster(clust_start, i, clust_num, version = version)
                ## reset cluster number
                if not same_chrom:
                    clust_num = 1
                else:
                    clust_num += 1
            if is_last:
                dat[i].append("singleton")
            clust_start = i
        elif (same_chrom and within_range) and is_last: ## same cluster
            assign_cluster(clust_start, i + 1, clust_num, version = version)
        last_entry = entry
    
    ## write to file
    output_colnames = ["chrom", "start", "end", "gene", "strand", "cluster"]
    f = open(fout(cluster_within_kb), "w+")
    f.write('\n'.join(['\t'.join(map(str, x))
                       for x in [output_colnames] + get(dat, *output_colnames)]) + '\n')
    f.close()
    return
