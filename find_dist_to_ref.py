import os
import re
import sys
sys.path.append("/mnt/chaelab/rachelle/src")
from basics import *
from data_manip import *
from Bio import Phylo

## presetss
REGEX_GID_ARATH = "AT.G\d+"
COI_F = "/mnt/chaelab/rachelle/hap_tags/results/nlr165/nlr165_final_wPred.tsv"

def write_closest_ref(t_nwk, ref_pattern, *ref_patterns,
                      fouts = {}, patterns = {}, ref_gid_pattern = "^.+$",
                      coi_f = None, coi_uncategorised = None,
                      outgroup_pattern = None, rooted = False) -> None:
    """
    Get and write closest reference gene for each leaf in tree.
    
    Arguments:
        t_nwk (str): required, path to nwk tree
        ref_pattern (str): required, regex pattern of leaf names of reference genes
        *ref_patterns (str): optional, additional regex pattern for further filtering after
            initial execution of <Bio.Phylo tree>.find_elements
        fouts (dict): required, dict of paths to output files, indexed identically to arg 'patterns',
            where data for leaf names matching the regex pattern of given index 
            will be written to the corresponding file
        patterns (dict): required, regex pattern of leaf names to assign, 
            indexed identically to arg 'fouts'
        ref_gid_pattern (str): optional, regex pattern to extract gene ID from reference gene leaf name
            (default='^.+$')
        coi_f (str): optional, path to cluster assignment file
        coi_uncategorised (str): optional, cluster name to assign to genes that are not covered by coi_f
            (default=None)
        outgroup_pattern (str): optional, regex pattern of root node,
            used if rooted=False
        rooted (bool): whether input tree is rooted (default=False);
            if False and outgroup_pattern=None, tree will be rooted arbitraily at midpoint
    """
    ## check that there's a 1-to-1 mapping of keys in fouts and patterns
    if (set(fouts.keys()) != set(patterns.keys())):
        print("1-to-1 mapping of keys in args 'fouts' and 'patterns' is required")
        return
    ## read tree
    t = Phylo.read(t_nwk, "newick")
    if not rooted:
        if not root_pattern:
            print("Rooting tree at midpoint")
            t.root_at_midpoint()
        elif not rooted:
            print("Rooting tree with outgroup")
            root = tuple(t.find_elements(name=root_pattern))
            t.root_with_outgroup(root)
    ## get reference nodes
    ref_nodes = tuple(t.find_elements(name=ref_pattern, terminal=True))
    for pattern in ref_patterns:
        ref_nodes = [n for n in ref_nodes if re.search(pattern, n.name)]
    collapse_iter = lambda iterable, char: [char.join([str(y) for y in x]) for x in iterable]
    def generate_to_write (nodes):
        to_write = '\n'.join(['\t'.join(["contig", "gene", "cluster",
                                         "distance", "ref_name", "sisters"])] +
                             ['\t'.join(collapse_iter([(node.name,), node.gene,
                                                       node.cluster, (node.closest_dist,),
                                                       tuple(r.name for r in node.closest_ref),
                                                       tuple(collapse_iter([(sis.name, d) for sis, d in
                                                                            node.ref_sisters.items()], ','))],
                                                      ';'))
                              for node in nodes])
        return to_write
    for k, pattern in patterns.items():
        print(f"Working on {k}")
        nodes = tuple(t.find_elements(name=pattern, terminal=True))
        assign_closest_gene_and_cluster(t, nodes, ref_nodes, gene_pattern = ref_gid_pattern,
                                        coi_f = coi_f, uncategorised = coi_uncategorised)
        with open(fouts[k], "w+") as f:
            f.write(generate_to_write(nodes))
    return

def assign_closest_gene_and_cluster(t, q_nodes, ref_nodes, gene_pattern = "^.+$",
                                    coi_f = None, uncategorised = "singletons"):
    if coi_f:
        clusters = read_clusters(coi_f)
    else:
        clusters = {}
    print("Assigning genes and clusters to reference nodes")
    assign_gene_cluster_to_ref(ref_nodes, clusters, gene_pattern = gene_pattern, uncategorised = uncategorised)
    print("Finding inner distance to references")
    assign_inner_dist_to_ref(t, ref_nodes)
    print("Finding distances from references to queries")
    assign_dist_to_ref(t, q_nodes)
    print("Finding closest reference to queries")
    for i, q_node in enumerate(q_nodes):
        closest_ref = min(q_node.ref_distance, key = lambda x: q_node.ref_distance[x])
        q_node.closest_dist = q_node.ref_distance[closest_ref]
        q_node.closest_ref = tuple(ref_node for ref_node, dist in q_node.ref_distance.items() if
                                   dist == q_node.closest_dist)
        q_node.gene = tuple(ref_node.gene for ref_node in q_node.closest_ref)
        q_node.cluster = tuple(ref_node.cluster for ref_node in q_node.closest_ref)
        # q_node.domain = domain
        # q_node.domain_name = tuple(ref_node.domain for ref_node in q_node.closest_ref)
        # q_node.isoform_domain = tuple(ref_node.isoform_domain for ref_node in q_node.closest_ref)
    return

def assign_dist_to_ref(t, q_nodes):
    root = t.root
    for i, q_node in enumerate(q_nodes):
        if i % 500 == 0:
            print("Query", i)
        parent_nodes = [root] + t.get_path(q_node)
        distance = q_node.branch_length
        q_node.ref_distance = {}
        q_node.ref_sisters = {}
        for parent_node in parent_nodes[-2::-1]:
            try:
                q_node.ref_distance = {**q_node.ref_distance,
                                       **{k: v + distance for k, v in parent_node.ref_distance.items()
                                          if not k in q_node.ref_distance}}
                if not q_node.ref_sisters:
                    q_node.ref_sisters = {k: v + distance for k, v in parent_node.ref_distance.items()}
            except AttributeError:
                continue
            try:
                distance += parent_node.branch_length
            except TypeError:
                continue
    return

def assign_inner_dist_to_ref(t, ref_nodes):
    root = t.root
    for ref_node in ref_nodes:
        parent_nodes = [root] + t.get_path(ref_node)
        distance = ref_node.branch_length
        for parent_node in parent_nodes[-2::-1]:
            try:
                parent_node.ref_distance = {**parent_node.ref_distance,
                                            **{ref_node: distance}}
            except AttributeError:
                parent_node.ref_distance = {ref_node: distance}
            try:
                distance += parent_node.branch_length
            except TypeError:
                continue
    return

def read_clusters(coi_f):
    ## with open("/mnt/chaelab/rachelle/hap_tags/data/nlr164_wPred.tsv", 'r') as f:
    ##    coi = [x[:-1].split('\t') for x in f.readlines()]
    # with open("/mnt/chaelab/rachelle/hap_tags/results/nlr165/nlr165_final_wPred.tsv", 'r') as f:
    #    coi = [x[:-1].split('\t') for x in f.readlines()]
    coi = [x.split('\t') for x in splitlines(coi_f)]
    get = make_custom_get(coi[0])
    coi = coi[1:]
    output = {cluster: {"genes": set()} for cluster in set(get(coi, "cluster"))}
    for entry in coi:
        output[get(entry, "cluster")]["genes"] |= {get(entry, "gene")}
    # with open(coi_f, 'r') as f:
    #     coi = [x[:-1].split(' ') for x in f.readlines()]
    # header = coi[0]
    # get = make_custom_get(header)
    # coi = coi[1:]
    # output = {}
    # for entry in coi:
    #     with open(os.path.join(cluster_dir, get(entry, "clust_name") + ".txt"), 'r') as f:
    #         output[get(entry, "dir_name")] = {**{colname: get(entry, colname) for colname in header},
    #                                             **{"genes": set(x[:-1] for x in f.readlines())}}
    return output

def assign_gene_cluster_to_ref(ref_nodes, clusters, gene_pattern = "^.+$", uncategorised = "singletons"):
    for ref_node in ref_nodes:
        ref_node.gene = re.search(gene_pattern, ref_node.name).group(0)
        # ref_node.domain = ref_node.gene + '.' + re.search("\d+(?=\|complete)", ref_node.name).group(0)
        # ref_node.isoform_domain = re.search("(?<=\|)AT.G\d+?.\d+?\|.+?\|\d+", ref_node.name).group(0)
        possible_clusters = tuple(k for k, v in clusters.items() if ref_node.gene in v["genes"])
        ref_node.cluster = uncategorised if not possible_clusters else possible_clusters[0]
    return

def get_closest_ref(dist_to_ref):
    output = {q_node: min(v, key = lambda x: v[x]) for q_node, v in dist_to_ref.items()}
    return output    

def get_dist_to_ref(t, q_nodes, ref_nodes):
    output = {q_node: {ref_node: t.distance(q_node, ref_node)
                       for ref_node in ref_nodes}
              for q_node in q_nodes}
    return output
