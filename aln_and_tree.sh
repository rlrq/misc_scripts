#!/bin/bash

## although mafft and FastTree are referred to as ALN and TREE, this script will ONLY work for mafft (cuz we're using the --add parameter to chaining alignments) and FastTree (cuz it needs < for input and not all tools support that) lol. ALN and TREE are mostly there in case there are different versions of these tools (e.g. FastTreeMP)
INPT=()
ALN_DEFAULT=/usr/bin/mafft
TREE_DEFAULT=/mnt/chaelab/rachelle/programmes/FastTree
PREFIX_DEFAULT=alnAndTree

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

## set some variables
DIR="${DIR:-$(pwd)}"
PREFIX="${PREFIX:-${PREFIX_DEFAULT}}"
ALN="${ALN:-${ALN_DEFAULT}}"
TREE="${TREE:-${TREE_DEFAULT}}"
if [ $(basename ${ALN}) == "mafft" ]; then
    aln_suff='mafft'
else
    aln_suff='aln'
fi
if [ $(basename ${TREE}) == "FastTree" ] || [ $(basename ${TREE}) == "FastTreeMP" ]; then
    tree_ext='nwk'
else
    tree_ext='tree'
fi
ALN_OUT="${ALN_OUT:-${DIR}/${PREFIX}.${aln_suff}.fasta}"
TREE_OUT="${TREE_OUT:-${DIR}/${PREFIX}.${aln_suff}.${tree_ext}}"

## print some stuff for checking
echo "Alignment output: ${ALN_OUT}"
echo "Tree output: ${TREE_OUT}"
for inpt in ${INPT[@]}; do
    echo ${inpt}
done

## argument checking
if [ ${#INPT[@]} -eq 0 ] && [ -z ${SEED} ]; then
    echo "At least one FASTA file is required (--seed <FASTA file of alignment> or positional arguments of unaligned FASTA file(s))"
    exit 1
fi

## if no seed file, align first file in INPT
if [ -z ${SEED} ]; then
    seed=$(mktemp -p ${DIR} "${PREFIX}.XXX" )
    ${ALN} ${ALN_ARGS} ${INPT[0]} > ${seed}
    INPT=( ${INPT[@]:1} ) ## remove first file from INPT
    ## move file to ALN_OUT (ALN_OUT will be overwritten if there are other fasta files to add)
    mv ${seed} ${ALN_OUT}
    SEED=${ALN_OUT}
    # ## if no other input files, move this alignment to ALN_OUT
    # if [ ${#INPT[@]} -eq 0 ]; then
    #     mv ${seed} ${ALN_OUT}
    #     SEED=${ALN_OUT}
    # else
    #     SEED=${seed}
    # fi
else
    ## if no other input files, treat SEED as ALN_OUT
    if [ ${#INPT[@]} -eq 0 ]; then
        ALN_OUT=${SEED}
    fi
fi

## while there are files to add to existing alignment, align
seed=${SEED}
for inpt in ${INPT[@]}; do
    echo "Aligning ${inpt}"
    tmp=$(mktemp -p ${DIR} "${PREFIX}.XXX")
    ${ALN} ${ALN_ARGS} --add ${inpt} ${seed} > ${tmp}
    mv ${tmp} ${ALN_OUT}
    seed=${ALN_OUT}
done

## proceed to tree building
${TREE} ${TREE_ARGS} < ${ALN_OUT} > ${TREE_OUT}

exit 0
