#!/usr/bin/bash

## uses process_OrthoFinder.py

## fout will be generated using summarise_mcl_output_v1 from test_programmes/src/summarise_mcl.py
rerun_OrthoFinder_mcl(){
    local dir prefix seqid_paired seqid_unpaired fout of_graph inflations mcl_args keep_tmp dir_tmp
    while [[ "${1}" ]]; do
        case "$1" in
            -d|--dir) dir=$(realpath ${2});; ## a temporary directory will be created within; temporary files will be moved to this directory if --keep-tmp is raised
            -p|--prefix) prefix="${2}";;
            -a|--algorithm) algorithm=${2};; ## blast/diamond/diamond_ultra_sens/blast_gz/mmseqs/blast_nucl
            --seqid-paired) seqid_paired=$(realpath ${2});; ## will be used to generate seqid_unpaired
            --seqid-unpaired) seqid_unpaired=$(realpath ${2});;
            -o|--out) fout=$(realpath ${2});;
            -g|--graph) of_graph=$(realpath ${2});; ## path to OrthoFinder_graph.txt in WorkingDirectory
            -I) inflations="${2}";; ## comma-separated inflation values to use
            --mcl) mcl=$(realpath ${2});; ## path to mcl
            --mcl-args) mcl_args="${2}";; ## miscellaneous mcl args
            --keep-tmp) keep_tmp=1;;
        esac
        shift
    done
    ## check args
    if [ -z "${algorithm}" ]; then
        echo "Algorithm needed. Use --algorithm <algo>. Algorithm options: blast, diamond, diamond_ultra_sens, blast_gz, mmseqs, blast_nucl"
        return 1
    fi
    ## set variables
    prefix=${prefix:-'OrthoFinder_regroup'}
    dir_tmp=$(mktemp -d)
    echo ${dir_tmp}
    ## make seqid_unpaired if not passed as input
    if [ -z "${seqid_unpaired}" ]; then
        if [ -z "${seqid_paired}" ]; then
            echo "SequenceIDs.txt needed if --seqid-unpaired is not used. Use --seqid-paired <path to SequenceIDs.txt>."
            return 1
        fi
        seqid_unpaired=${dir_tmp}/SequenceIDs.unpaired.txt
        python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src/'); import process_OrthoFinder; process_OrthoFinder.make_unpaired_id_map('${seqid_paired}', '${seqid_unpaired}')"
    fi
    ## run mcl for each inflation value
    dir_mcl=${dir_tmp}/mcl
    mkdir ${dir_mcl}
    inflation_array=$(echo ${inflations} | tr ',' ' ')
    for inflation in ${inflation_array[@]}; do
        ${mcl} ${of_graph} -I ${inflation} -o ${dir_mcl}/${prefix}.OrthoFinder-${algorithm}.mcl.${inflation}.txt -use-tab ${seqid_unpaired} ${mcl_args}
    done
    ## summarise
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/test_programmes/src/'); import summarise_mcl; summarise_mcl.summarise_mcl_output_v1('${prefix}', '${dir_mcl}', '${fout}')"
    ## cleanup
    if ! [ -z ${keep_tmp} ] && ! [ -z "${dir}" ] ; then
        echo "Moving files"
        mv ${dir_tmp}/* ${dir}
        rmdir ${dir_tmp}
    else
        echo "Removing temporary files"
        rm -r ${dir_tmp}
    fi
}

