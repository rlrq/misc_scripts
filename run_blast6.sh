#!/bin/bash

run_blast6 () {
    params=${@}
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            -outfmt) local executable="$(echo ${params} | sed "s/${2}/\x27${2}\x27/")";
                     local header="$(echo ${2} | sed 's/^[[:digit:]]\+ //' | tr ' ' '\t')";;
            -out) local fout="${2}";
                  local tmp="${fout}_tmp.txt";;
        esac
        shift
    done
    eval "${executable}" ## execute
    printf "${header}\n$(grep -v '^#' ${fout})" > $tmp; mv $tmp $fout
}

## run blastn with default parameters
## arg1: database
## arg2: output directory
## arg3: output prefix
## arg4: input fasta file of subject sequences
## arg5: output file name (optional)
blastn_default () {
    local db=${1}
    local dir=${2}
    local prefix=${3}
    local fasta=${4}
    local out_root=${dir}/${prefix}
    local fout=${5:-${out_root}.blastn_summary.tsv}
    local tmp=${out_root}_tmp.txt
    blastn -outfmt 11 -db ${db} -query ${fasta} -out ${fout} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe"
    printf "$(echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe' | tr ' ' '\t')\n$(grep -v '^#' ${fout})" > $tmp; mv $tmp $fout
}

## run discontiguous megablast
blastn_dc () {
    local db=${1}
    local dir=${2}
    local prefix=${3}
    local fasta=${4}
    local out_root=${dir}/${prefix}
    local fout=${5:-${out_root}.blastn_summary.tsv}
    local tmp=${out_root}_tmp.txt
    blastn -outfmt 11 -db ${db} -query ${fasta} -out ${results_asn} -task dc-megablast -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe"
    printf "$(echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe' | tr ' ' '\t')\n$(grep -v '^#' ${fout})" > $tmp; mv $tmp $fout
}

## run blastn instead of megablast
blastn_blastn () {
    local db=${1}
    local dir=${2}
    local prefix=${3}
    local fasta=${4}
    local out_root=${dir}/${prefix}
    local fout=${5:-${out_root}.blastn_summary.tsv}
    local tmp=${out_root}_tmp.txt
    blastn -outfmt 11 -db ${db} -query ${fasta} -out ${results_asn} -task blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe"
    printf "$(echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe' | tr ' ' '\t')\n$(grep -v '^#' ${fout})" > $tmp; mv $tmp $fout
}

## blast against Anna-Lena's 70 accessions Ren-Seq database
## arg1: output dir
## arg2: output prefix
## arg3: fasta file of subjects
## arg4: output file name (optional)
blastn_to_al70 () {
    blastn_default '/mnt/chaelab/shared/blastdb/anna_lena70.contigs/anna_lena70.contigs' ${1} ${2}.annaLena70 ${3} ${4}
}

## blast against known Col-0 NLRs
## arg1: output directory
## arg2: output prefix
## arg3: input fasta file of subject sequences
## arg4: output file name (optional)
blastn_to_col0nlr() {
    blastn_default '/mnt/chaelab/shared/blastdb/nlr164/all_NLR.164.db' ${1} ${2}.col0nlr ${3} ${4}
}
