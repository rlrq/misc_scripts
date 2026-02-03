#!/bin/bash

while (( "$#" )); do
    case "$1" in
        --gff) GFF="${2}"; shift;; ## input GFF3 file
        --blast-args) BLAST_ARGS="${2}"; shift;; ## arguments for blast
        --fasta) FASTA="${2}"; shift;; ## input FASTA file
        --id) F_FID="${2}"; shift;; ## input TXT file of feature IDs (1 ID per line)
        --extend) EXTEND_BP="${2}"; shift;; ## number of bp to extend for feature subsetting
        --flank) FLANK_BP="${2}"; shift;; ## number of bp to extend beyond extended features for FASTA subsetting
        --gff-out) GFF_OUT="${2}"; shift;; ## output GFF3 file
        --fasta-out) FASTA_OUT="${2}"; shift;; ## output FASTA file
        ## output BED file containing mapping coords to sequence IDs in FASTA_OUT
        --bed-out) BED_OUT="${2}"; shift;;
        ## write GFF3 subset with original coords to this file
        --gff-original|--gff-og) GFF_ORIGINAL="${2}"; shift;;
        ## move temporary directory to DIR_TMP_OUT upon exit instead of deleting it
        --tmp-out) DIR_TMP_OUT="${2}"; shift;;
        ## extend around coordinates of highest parent for each feature, not the feature itself
        --parent) EXTEND_PARENT=1;;
        --self) EXTEND_PARENT=0;; ## extend around coordinates of feature, not its highest parent
        --memsave) MEMSAVE=1;; ## use memsave mode
        --no-memsave) MEMSAVE=0;; ## don't use memsave mode
        ## sorts GFF3 file and then pass -sorted argument and sorted GFF3 file to bedtools intersect
        --sort) SORT=1;;
        ## chromosome size threshold for automatically setting SORT=1
        --sort-threshold|--threshold) THRESHOLD="${2}"; shift;;
    esac
    shift
done

