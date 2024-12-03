#!/bin/bash

## note: this uses Cdd.v3.18

while (( "$#" )); do
    case "$1" in
        -d|--dir) DIR="${2}"; shift;; ## output directory
        --pref|--prefix) PREFIX="${2}"; shift;; ## output prefix, used if --aln-out and/or --tree-out are not specified
        --pep|--fasta) FA_PEP="${2}"; shift;; ## fasta file of peptides to search for nlr domains for
    esac
    shift
done

DIR=${DIR:-$(pwd)}
PREFIX=${PREFIX:-rpsNLRdomain}

if [[ -z ${FA_PEP} ]]; then
    echo "Fasta file of peptide sequences required. Please provide using '--fasta <file>' or '--pep <file>'."
    exit 1
fi

DIR_SRC=/mnt/chaelab/rachelle/src

# ARATH_PEP_DB=/mnt/chaelab/shared/blastdb/tair10pep/tair10pep
# ARATH_NLR_PEP_FA=/mnt/chaelab/rachelle/data/NLR/nlr166.pep.fasta
# ARATH_NLRR_PEP_FA=/mnt/chaelab/rachelle/data/NLR/nlr166.rpw8.pep.fasta ## with RPW8 genes
# GID_FRM_ARATH_PEP_DB_PY="lambda x: '.'.join(x.split('.')[:-1])"
# GID_FRM_ARATH_NLR_PEP_FA_PY="lambda x: x.split('|')[1]"
# DIR_BRAPA_DAT=/mnt/chaelab/rachelle/data/Brapa/v3.0
# BRAPA_PEP_FA=${DIR_BRAPA_DAT}/Brapa_genome_v3.0_pep.fasta

# PEP_FA=${ARATH_NLR_PEP_FA}
# PEP_FA=${ARATH_NLRR_PEP_FA}
# PEP_DB=${ARATH_PEP_DB}
# GID_FRM_PEP_FA_PY=${GID_FRM_ARATH_NLR_PEP_FA_PY} ## targets (i.e. recip blast best hit must be to one of these)
# GID_FRM_PEP_DB_PY=${GID_FRM_ARATH_PEP_DB_PY} ## targets + background (for recip blast)
# QUERY_NAME=ArathNLR ## used by find_nlrs_in_other_indv to name output file of first_blast in format <prefix>.q<QUERY_NAME>.tsv

RPSBLAST=/usr/bin/rpsblast+
DB_CDD=/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd
CDD_VERSIONS=/mnt/chaelab/rachelle/data/cdd/2020_10_22/cdd.versions
declare -A DOMAINS
DOMAINS["TIR"]='366714';DOMAINS["RX-CC_like"]='271353';DOMAINS["RPW8"]='368547';DOMAINS["NB-ARC"]='366375';DOMAINS["Rx_N"]='375519';DOMAINS["LRR"]='227223';DOMAINS["PLN03210"]='215633'

# ## db required for num_threads in blastp
# build_tmp_db() {
#     local fasta=${1}
#     local iout=${2}
#     # echo "Building subject database"
#     makeblastdb -in ${fasta} -dbtype prot -out ${iout} -hash_index
# }

# ## requires FASTA file of all protein sequences in target indv
# first_blast() {
#     local sdb=${1}
#     local fout=${2}
#     local thread=${3:-4} ## default 4 threads
#     # echo "BLASTP A. thaliana NLR to subject"
#     blastp -query ${PEP_FA} -db ${sdb} -out ${fout} -outfmt 6 -num_threads ${thread}
# }

# get_candidate_seq() {
#     local blastout=${1} ## blast output of first_blast
#     local fasta=${2}
#     local fout=${3}
#     local fout_cds=${4}
#     # echo "Extracting candidate sequences"
#     local tmp_candidates_txt=$(mktemp)
#     echo "tmp file: ${tmp_candidates_txt}"
#     cut -f2 ${blastout} | sort | uniq > ${tmp_candidates_txt}
#     python3 -c "import sys; from Bio import SeqIO; sys.path.append('${DIR_SRC}'); from data_manip import splitlines; candidates = set(splitlines('${tmp_candidates_txt}')); to_write = {record.id: record.seq for record in SeqIO.parse('${fasta}', 'fasta') if record.id in candidates}; from fasta_manip import dict_to_fasta; dict_to_fasta(to_write, '${fout}')"
#     rm ${tmp_candidates_txt}
# }

# reciprocal_blast() {
#     local qcandidates=${1}
#     local fout=${2}
#     local threads=${3:-4}
#     # echo "Reciprocal BLASTP of candidate NLRs to A. thaliana"
#     blastp -query ${qcandidates} -db ${PEP_DB} -out ${fout} -outfmt 6 -num_threads ${thread}
# }

# filter_blast() {
#     local blastout=${1}
#     local fout=${2}
#     local topn=${3:-5}
#     python3 -c "import sys; sys.path.append('${DIR_SRC}'); from top_recipblast import blast6_top_n_memsave; blast6_top_n_memsave('${blastout}', fout = '${fout}', n = ${topn})"
# }

# get_nlrs_from_recip() {
#     local rblastout=${1}
#     local fout=${2}
#     ## comparepy is an anonymous python function to compare whether ARATH_PEP_DB's pep gene (db_id) == NLR_PEP_FA's pep gene (q_id)
#     local comparepy=${3:-"lambda q, s: (${GID_FRM_PEP_DB_PY})(q) == s"}
#     # echo "Getting IDs of potential NLRs"
#     python3 -c "import sys; sys.path.append('${DIR_SRC}'); from top_recipblast import filter_from_blast6; from fasta_manip import fasta_to_dict; targets = set(map(${GID_FRM_PEP_FA_PY}, fasta_to_dict('${PEP_FA}').keys())); filter_from_blast6('${rblastout}', '${fout}', targets = targets, eq = ${comparepy})"
# }

# filter_fasta() {
#     local gids=${1}
#     local fasta=${2}
#     local fout=${3}
#     python3 -c "import sys; sys.path.append('${DIR_SRC}'); from data_manip import splitlines; from fasta_manip import fasta_to_dict, dict_to_fasta; gids = set(splitlines('${gids}')); to_write = {k: v for k, v in fasta_to_dict('${fasta}').items() if k in gids}; dict_to_fasta(to_write, '${fout}')"
# }

nlr_rpsblast() {
    local query=${1}
    local fout=${2}
    ${RPSBLAST} -query ${query} -db ${DB_CDD} -out ${fout} -outfmt 6
}

# ## e.g. ~$ aarray_to_dictpy DOMAINS
# aarray_to_dictpy() {
#     ## convert associative array to python dictionary format (all strings)
#     local -n array=${1} ## just the array name
#     local output=''
#     for k in ${!array[@]}; do
#         local v=${array[${k}]}
#         local output="${output} '${k}': '${v}',"
#     done
#     echo "{${output}}"
# }

filter_rpsblast() {
    local rpsblastout=${1}
    local fout=${2}
    python3 -c "import sys; sys.path.append('${DIR_SRC}'); from rpsblast_manip import filter_rpsblast; filter_rpsblast('${rpsblastout}', '${fout}', '$(echo ${DOMAINS[@]})'.split(' '), cols_to_write = (0, 6, 7, 1), col_names_out = ('qseqid', 'qstart', 'qend', 'sseqid', 'domain'), write_domain_name = True, cdd_versions = '${CDD_VERSIONS}')"
}

get_nlr_architecture() {
    local rblastfilt=${1} ## output of get_nlrs_from_recip
    local fout=${2}
    local code="{'T': (366714,), 'C': (271353, 375519), 'R': (368547,), 'N': (366375,), 'L': (227223,)}"
    local order="['C', 'T', 'R', 'N', 'L']"
    python3 -c "import sys; sys.path.append('${DIR_SRC}'); from rpsblast_manip import rps_to_nlr; rps_to_nlr('${rblastfilt}', '${fout}', ${code}, by_start = True, collapse = True, collapse_order = ${order}, merge_overlapping_with_same_code = True)"
}

## function to execute RPS-BLAST and get NLR domain architecture
rpsblast_and_get_nlr_architecture() {
    local outdir=${1}
    local prefix=${2}
    local pep=${3}
    mkdir -p ${outdir}
    echo "RPSBLAST"
    local filtrps=${outdir}/${prefix}.rpsblast.tsv
    nlr_rpsblast ${pep} ${filtrps}
    echo "Getting candidates with NB-ARC, TIR, RPW8, RX-CC_like, Rx_N, or LRR domains"
    local nlr_domains=${outdir}/${prefix}.NLR.domains.tsv
    filter_rpsblast ${filtrps} ${nlr_domains}
    echo "Getting NLR gids"
    local nlr_gid=${outdir}/${prefix}.NLR.gid.txt
    cut -f1 ${nlr_domains} | awk 'NR>1' | sort | uniq > ${nlr_gid}
    echo "Summarising architecture of NLR by domains and closest A. thaliana homologue"
    local nlr_architecture=${outdir}/${prefix}.NLR.architecture.tsv
    get_nlr_architecture ${filtrps} ${nlr_architecture}
}

rpsblast_and_get_nlr_architecture ${DIR} ${PREFIX} ${FA_PEP}
