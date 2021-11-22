from Bio import Phylo

## wrapper for collapse_clade_by_branch_length
def collapse_tree_by_branch_length_from_file(tree, tree_type, branch_threshold = 0):
    t = Phylo.read(tree, tree_type)
    collapse_tree_by_branch_length(t, branch_threshold)
    return t

def collapse_tree_by_branch_length(t, branch_threshold = 0):
    r = t.root
    collapse_clade_by_branch_length(r, branch_threshold)
    return None

## ## if (1) node A has children nodes B and C, (2) node B has children nodes D and E, and
## (3) the length of the branch between A and B <= branch_threshold,
## node B will be removed from node A's children and children nodes of B (D and E)
## will be added to node A's direct children, such that the node A will have children nodes C, D, and E.
def collapse_clade_by_branch_length(parent_clade, branch_threshold = 0):
    for clade in parent_clade:
        if not clade.is_terminal():
            ## collapse any sub-clades before current clade to avoid accidentally skipping any branches
            collapse_clade_by_branch_length(clade, branch_threshold)
            if clade.branch_length <= branch_threshold:
                parent_clade.collapse(clade)
    return None

## work on all (newick) trees in a directory
def collapse_trees_in_dir_by_branch_length_nwk(directory, tree_type = "newick", branch_threshold = 0,
                                               out_dir = "", in_place = False, act_on_subdirs = False):
    if not in_place:
        if not out_dir:
            print("Please provide an output directory if not executing this function in-place.")
            return 1
        elif not os.path.exists(out_dir):
            os.makedirs(out_dir)
    for fname in os.listdir(directory):
        if tree_type == "newick" and fname.split('.')[-1] == "nwk":
            t = collapse_tree_by_branch_length(os.path.join(directory, fname), tree_type, branch_threshold)
            Phylo.write(t, os.path.join((directory if in_place else out_dir), fname), tree_type)
    # for subdirs, dirs, files in os.walk(directory):
    #     for file in files:
    #         t = collapse_tree_by_branch_length(os.path.join(subdir[0], file), tree_type, branch_threshold)
    #         Phylo.write(t, os.path.join(directory if in_place else out_dir, file) + "newick")
    #     if act_on_subdirs:
    #         for i, subdir in enumerate(subdirs[1:]):
    #             collapse_trees_in_dir_by_branch_length(subdir, tree_type, branch_threshold, out_dir = subdir if in_place else os.path.join(out_dir, dirs[i + 1]), in_place = in_place, act_on_subdirs = act_on_subdirs)
    return 0
