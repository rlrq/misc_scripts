#!/bin/bash

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
    local results_asn=${out_root}.blastn.asn
    local results_sum=${out_root}.blastn_summary.txt
    blastn -outfmt 11 -db ${db} -query ${fasta} -out ${results_asn}
    blast_formatter -archive ${results_asn} -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe" > ${results_sum}
    printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\ts. len\tevalue\tbit score\tq. frame\ts. frame\n$(grep -v '^#' ${results_sum})" > $fout
}

## run discontiguous megablast
blastn_dc () {
    local db=${1}
    local dir=${2}
    local prefix=${3}
    local fasta=${4}
    local out_root=${dir}/${prefix}
    local fout=${5:-${out_root}.blastn_summary.tsv}
    local results_asn=${out_root}.blastn.asn
    local results_sum=${out_root}.blastn_summary.txt
    blastn -outfmt 11 -db ${db} -query ${fasta} -out ${results_asn} -task dc-megablast
    blast_formatter -archive ${results_asn} -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe" > ${results_sum}
    printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\ts. len\tevalue\tbit score\tq. frame\ts. frame\n$(grep -v '^#' ${results_sum})" > $fout
}

## run blastn instead of megablast
blastn_blastn () {
    local db=${1}
    local dir=${2}
    local prefix=${3}
    local fasta=${4}
    local out_root=${dir}/${prefix}
    local fout=${5:-${out_root}.blastn_summary.tsv}
    local results_asn=${out_root}.blastn.asn
    local results_sum=${out_root}.blastn_summary.txt
    blastn -outfmt 11 -db ${db} -query ${fasta} -out ${results_asn} -task blastn
    blast_formatter -archive ${results_asn} -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore qframe sframe" > ${results_sum}
    printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\ts. len\tevalue\tbit score\tq. frame\ts. frame\n$(grep -v '^#' ${results_sum})" > $fout
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
