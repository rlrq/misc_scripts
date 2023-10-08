#!/bin/bash

DIR_SRC=/mnt/chaelab/rachelle/src
FastTree=/mnt/chaelab/rachelle/programmes/FastTree
newick=/mnt/chaelab/rachelle/programmes/newick/src/o/newick
newick_tools=/mnt/chaelab/rachelle/programmes/newick-tools/src/newick-tools

source ${DIR_SRC}/getseq_manip.sh

get_homologue_genomic() {
    local fout=${1}
    local gid=${2}
    local directories=${@:3}
    ## clear output file
    touch ${fout}
    rm ${fout}
    ## tempfile
    local tmp=$(mktemp)
    touch ${tmp}
    echo ${tmp}
    for dir in ${directories[@]}; do
        ## for each of the relevant file in the directory (should only be one, but we're playing safe)
        for f in ${dir}/final/*.final.gene.fasta; do
            ## get homologue IDs
            local txt_report=${f%*.final.gene.fasta}.report.txt
            local homologues=$(awk -v gid=${gid} '$2==gid' ${txt_report} | cut -f1 | tr '\n' ',')
            ## get homologue sequences
            get_seq_from_gid ${f} ${tmp} ${homologues} gene
            ## cat to final
            cat ${tmp} >> ${fout}
        done
    done
    rm ${tmp}
}

get_closest_nonoptimal_homologue_genomic() {
    local fout=${1}
    local gid=${2}
    local directories=${@:3}
    ## clear output file
    touch ${fout}
    rm ${fout}
    ## tempfile
    local tmp=$(mktemp)
    touch ${tmp}
    echo ${tmp}
    local dir
    for dir in ${directories[@]}; do
        ## for each of the relevant file in the directory (should only be one, but we're playing safe)
        local f
        for f in ${dir}/final/*.final.gene.fasta; do
            ## get non-optimal homologue entries, append to tmp
            local txt_report=${f%*.final.gene.fasta}.report.txt
            awk -v gid=${gid} -v d=${dir} -v fname=${f} -v OFS='\t' '$4!=1 && $5==gid {print $0,d,fname}' ${txt_report} >> ${tmp}
        done
    done
    ## get best non-optimal homologue
    local entry=$(sort -k6n ${tmp} | tail -1)
    local homologue=$(printf "${entry}" | cut -f1)
    local homologue_dir=$(printf "${entry}" | cut -f7)
    local homologue_fname=$(printf "${entry}" | cut -f8)
    local homologue_pref=$(basename ${homologue_fname%*.final.gene.fasta})
    get_seq_from_gid ${homologue_dir}/blast1/${homologue_pref}.blast1.*.intersect_query.gene.fasta \
                     ${fout} ${homologue} gene
    rm ${tmp}
}

align_and_tree_homologue_and_nonoptimal_genomic() {
    ## outgrp should contain EXACTLY ONE record
    local fout_aln=${1}
    local fout_nwk=${2}
    local gid=${3}
    local outgrp=${4}
    local directories=${@:5}
    ## clear output file
    touch ${fout_aln} ${fout_nwk}
    rm ${fout_aln} ${fout_nwk}
    ## get homologues
    local fa_optimal=$(mktemp)
    local fa_nonoptimal=$(mktemp)
    get_homologue_genomic ${fa_optimal} ${gid} ${directories[@]}
    get_closest_nonoptimal_homologue_genomic ${fa_nonoptimal} ${gid} ${directories[@]}
    ## if no sequence in ${nonoptimal_fa}, use the outgroup file provided by user
    if [ $(grep '>' ${fa_nonoptimal} | wc -l) -eq 0 ]; then
        cp ${outgrp} ${fa_nonoptimal}
    fi
    ## align
    local fa_aln=$(mktemp)
    local nwk_tree=$(mktemp)
    mafft ${fa_optimal} > ${fa_aln}
    mafft --add ${fa_nonoptimal} ${fa_aln} > ${fout_aln}
    ## tree
    ${FastTree} -nt ${fout_aln} > ${fout_nwk}
    ## root by outgroup
    local outgrp_id="$(grep -Po '(?<=^>).+' ${fa_nonoptimal})"
    echo ${outgrp_id}
    ${newick_tools} --tree_file ${fout_nwk} --root "${outgrp_id}" --output_file ${nwk_tree}
    ## remove outgroup
    ${newick_tools} --tree_file ${nwk_tree} --prune_tips "${outgrp_id}" --output_file ${fout_nwk}
    ## remove tmp files
    rm ${fa_optimal} ${fa_nonoptimal} ${fa_aln} ${nwk_tree}
}

rename_tree() {
    ## outgrp should contain EXACTLY ONE record
    local fin=${1}
    local fout=${2}
    local rename=${3}
    local outgrp=${4}
    local tmp=$(mktemp)
    cp ${rename} ${tmp}
    if [ ! -z "${outgrp}" ]; then
        ## append outgroup's new label to a copy of $rename
        local outgrp_id="$(grep -Po '(?<=^>).+' ${outgrp})"
        printf "${outgrp_id}\toutgroup.outgroup\n" >> ${tmp}
    fi
    ${DIR_SRC}/rename_newick.py ${fin} ${tmp} ${fout}
    ## remove $rename copy
    rm ${tmp}
}

collapse_tree_by_phyla() {
    local fin=${1}
    local fout=${2}
    local features=${3}
    ${newick} -condense ${fin} -features ${features} -output ${fout}
}

rofo_of_selected_phyla() {
    
}
