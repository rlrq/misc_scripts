#!/bin/bash

## if INPT is empty and SEED is provided, ALN_OUT will not be written (SEED will be treated as ALN_OUT)
## shift twice after a keyword argument to skip over its value because otherwise anything that isn't a keyword argument will be added to INPT
while (( "$#" )); do
    case "$1" in
        --dir) DIR="${2}"; shift;; ## output directory
        --pref|--prefix) PREFIX="${2}"; shift;; ## output prefix, used if --aln-out and/or --tree-out are not specified
        --seed) SEED="${2}"; shift;; ## aligned seqs to seed further alignment instead of building alignment from scratch
        --aln) ALN_ARGS="${2}"; shift;; ## arguments for alignment
        --mafft) ALN_ARGS="${2}"; shift;; ## arguments for alignment
        --tree) TREE_ARGS="${2}"; shift;; ## arguments for tree building
        --fasttree) TREE_ARGS="${2}"; shift;; ## arguments for tree building
        --aln-out) ALN_OUT="${2}"; shift;; ## output alignment file
        --tree-out) TREE_OUT="${2}"; shift;; ## output tree file
        --which-aln) ALN="${2}"; shift;; ## path to alignment binary
        --which-tree) TREE="${2}"; shift;; ## path to tree building binary
        *) INPT+=( "${1}" ); ## input files will be progressively added to alignment order
    esac
    shift
done

## yknow what just use aln_and_tree with like extra options
