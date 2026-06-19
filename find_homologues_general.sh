#!/bin/bash

## based on find_homologues_in_phytozome.sh
## but generalised to take a lookup file such as in /mnt/chaelab/rachelle/data/genome_sets/
## to define genome fa, gff, and pep fa files instead of downloading from phytozome (requires header line "##@fields:...")

## Takes genome/proteome ID(s) and lookup file
## Takes Arabidopsis AGI gene ID(s)
## Extracts query Arabidopsis protein sequence(s)
## Reads location of genome's proteome, genome, and GFF file from lookup file
## Executes reciprocal protein blast on proteome using Arabidopsis protein as query
## Maps homologue transcript ID(s) to gene ID(s)
## Extracts GFF of relevant gene(s) from genomes
## Subsets genome according to GFF
## Optionally generates genbank file(s)

rblastgen=/mnt/chaelab/rachelle/scripts/recip_blast/recip_blast_gen.sh
subset_gff=/mnt/chaelab/rachelle/src/subset_gff.sh
get_seq=/mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh

GFF_DEFAULT=/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.gff
FASTA_DEFAULT=/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta
PEP_DEFAULT=/mnt/chaelab/shared/genomes/TAIR10/fasta/TAIR10_pep_20101214.fasta
# GENOME_LOOKUP_DEFAULT=/mnt/chaelab/rachelle/data/Phytozome/phytozome14/Current_Genomes.unrestricted.assembly-gff-prot_available-standardised.curated.tsv
PROTEOME_ID_DEFAULT=all
PROCESSES_DEFAULT=1
EXTEND_BP_DEFAULT=0
FLANK_BP_DEFAULT=500
GENBANK_OUT_DEFAULT=0
# EXTEND_PARENT_DEFAULT=1
MEMSAVE_DEFAULT=1
SORT_DEFAULT=0
THRESHOLD_DEFAULT=511999999 ## maximum chromosome size for bedtools without -sorted option is 512 Mb (https://bedtools.readthedocs.io/en/latest/)

while (( "$#" )); do
    case "$1" in
        --gff) GFF="$(realpath ${2})"; shift;; ## input GFF3 file of reference organism (TAIR10 default)
        --fasta) FASTA="$(realpath ${2})"; shift;; ## input FASTA file of reference assembly (TAIR10 default)
        --pep) PEP="$(realpath ${2})"; shift;; ## input FASTA file of reference proteome (TAIR10 default)
        --query-pep) PEP_QUERY="$(realpath ${2})"; shift;; ## input FASTA file of query protein sequences (not compatible with GID and F_GID)
        --id-file) F_GID="$(realpath ${2})"; shift;; ## input TXT file of gene IDs (1 ID per line)
        --id) GID="${2}"; shift;; ## comma-separated gene IDs
        --genome-lookup) GENOME_LOOKUP="$(realpath ${2})"; shift;; ## file with genome IDs in 1st column
        -p|--proteome) PROTEOME_ID="${2}"; shift;; ## comma-separated genome/proteome IDs
        --process) PROCESSES="${2}"; shift;; ## number of parallel processes
        --extend) EXTEND_BP="${2}"; shift;; ## number of bp to extend for feature subsetting
        --flank) FLANK_BP="${2}"; shift;; ## number of bp to extend beyond extended features for FASTA subsetting
        --dir) DIR="$(realpath ${2})"; shift;; ## output directory
        --gb|--genbank) GENBANK_OUT=1;; ## output genbank files
        --memsave) MEMSAVE=1;; ## use memsave mode
        --no-memsave) MEMSAVE=0;; ## don't use memsave mode
        ## sorts GFF3 file and then pass -sorted argument and sorted GFF3 file to bedtools intersect
        --sort) SORT=1;;
        ## chromosome size threshold for automatically setting SORT=1
        --sort-threshold|--threshold) THRESHOLD="${2}"; shift;;
    esac
    shift
done

## set some variables with defaults
GFF="${GFF:-${GFF_DEFAULT}}"
FASTA="${FASTA:-${FASTA_DEFAULT}}"
PEP="${PEP:-${PEP_DEFAULT}}"
GENOME_LOOKUP="${GENOME_LOOKUP:-${GENOME_LOOKUP_DEFAULT}}"
PROTEOME_ID="${PROTEOME_ID:-${PROTEOME_ID_DEFAULT}}"
GENBANK_OUT="${GENBANK_OUT:-${GENBANK_OUT_DEFAULT}}"
PROCESSES="${PROCESSES:-${PROCESSES_DEFAULT}}"
FLANK_BP="${FLANK_BP:-${FLANK_BP_DEFAULT}}"
EXTEND_BP="${EXTEND_BP:-${EXTEND_BP_DEFAULT}}"
MEMSAVE=${MEMSAVE:-${MEMSAVE_DEFAULT}}
SORT=${SORT:-${SORT_DEFAULT}}
THRESHOLD=${THRESHOLD:-${THRESHOLD_DEFAULT}}

## check that required arguments are provided
if [ -z ${GFF} ]; then
    echo "Input GFF3 file is required."
    exit 1
elif [ -z ${FASTA} ]; then
    echo "Input FASTA file is required."
    exit 1
elif [ -z ${PEP_QUERY} ] && [ -z ${F_GID} ] && [ -z ${GID} ]; then
    echo "Input protein FASTA file or gene ID file (1 ID per line) or string of comma-separated gene IDs is required."
    exit 1
elif [ ! -z ${PEP_QUERY} ] && ( [ ! -z ${F_GID} ] || [ ! -z ${GID} ] ); then
    echo "Input protein FASTA file (--pep-query) is incompatible with gene ID file (--id-file) and string of gene IDs (--id)."
    exit 1
elif [ -z ${GENOME_LOOKUP} ]; then
    echo "Genome lookup file is required."
    exit 1
elif [ -z ${DIR} ]; then
    echo "Output directory is required."
    exit 1
fi

## create temporary directory
mkdir -p ${DIR}
dir_tmp=${DIR}/tmp
mkdir -p ${dir_tmp}
echo "============::=== Working directory: '${DIR}' ==::============"
function mktemp_f { echo "$(mktemp -p ${dir_tmp})"; }

## cleanup
# function cleanup_delete { echo "Removing temporary directory '${dir_tmp}'"; rm -r ${dir_tmp}; }
# function cleanup_move { echo "Moving temporary directory '${dir_tmp}' to '${DIR}'";
#                         rm -r ${dir_tmp}/tmp; mv ${dir_tmp} ${DIR}; }
# if [ -z ${DIR} ]; then ## since we require DIR, we shouldn't encounter this
#     trap cleanup_delete EXIT
#     trap cleanup_delete SIGINT
# else
#     trap cleanup_move EXIT
#     trap cleanup_move SIGINT
# fi
## since we're writing directly to output directory, just delete the tmp directory inside it
function cleanup_delete_tmp { echo "Removing temporary directory '${dir_tmp}'"; rm -r ${dir_tmp}; }
trap cleanup_delete_tmp EXIT
trap cleanup_delete_tmp SIGINT

## get query protein sequences
if [ -z ${PEP_QUERY} ]; then
    ## create a copy of F_GID in our working directory
    mkdir -p ${DIR}/id
    if [ ! -z ${F_GID} ]; then
        cp ${F_GID} ${DIR}/id/query.id.txt
    fi
    if [ ! -z ${GID} ]; then
        if [ ! -z ${F_GID} ]; then
            echo "Both --id-file and --id were provided. Combined IDs will be used."
        else
            printf '' > ${DIR}/id/query.id.txt
        fi
        echo "${GID}" | tr ',' '\n' >> ${DIR}/id/query.id.txt
        F_GID=${DIR}/id/query.id.txt
    fi
    ## extract query sequences
    echo "Extracting query sequences"
    mkdir -p ${DIR}/seq
    ${get_seq} --assembly ${FASTA} --gff ${GFF} --feature CDS \
               --gene "$(cat "${F_GID}" | tr '\n' ',' | sed 's/,\+$//')" \
               --adjust-dir --translate --out ${DIR}/seq/query.pep.fasta
else
    mkdir -p ${DIR}/seq
    cp ${PEP_QUERY} ${DIR}/seq/query.pep.fasta
fi
PEP_QUERY=${DIR}/seq/query.pep.fasta

get_genome_lookup_v1 () {
    local proteome_id=${1}
    local field=${2}
    awk -v pid=${proteome_id} '( $1==pid || $1~/##@fields:/)' ${GENOME_LOOKUP} |
        sed 's/##@fields://' |
        awk -v col=${field} -v FS='\t' 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR>1{print $c}'
}

generate_fids_v1 () {
    echo "function: generate_fids_v1"
    local f_id=${1}
    local fout=${2}
    local patterns="${3}"
    python3 -c "
import itertools, re
patterns = set(re.search('(?<==)[^=]+', x).group(0) for x in '${patterns}'.split(';'))
print(patterns)
pattern_substitutions = {pattern: re.findall('(?<=<)[^>]+', pattern) for pattern in patterns}
substitutions = set(itertools.chain(*pattern_substitutions.values()))
grep_patterns = {s: ('^[^ ]+' if s == '0' else ('(?<=' + s + '=)[^ ]+')) for s in substitutions}
print(grep_patterns)
with open('${f_id}', 'r') as fin:
    with open('${fout}', 'w+') as fout:
        for seqid in fin:
            seqid = seqid.rstrip()
            values = {s: re.search(grep_pattern, seqid).group(0) for s, grep_pattern in grep_patterns.items()}
            for pattern in patterns:
                p = pattern
                for pattern_sub in pattern_substitutions[pattern]:
                    p = p.replace('<' + pattern_sub + '>', values[pattern_sub])
                if pattern != '':
                    _ = fout.write(p + '\n')
"
}

## find potential homologues
find_general_homologues_v1 () {
    echo "function: find_general_homologues_v1"
    local proteome_id=${1}
    local dir_out=${2}
    local fa_pep_query=${3}
    local fa_pep_bg=${4}
    local fa_pep=${5}
    local fa_assembly=${6}
    local fa_gff=${7}
    local map_patterns="${8}"
    mkdir -p ${dir_out} ${dir_out}/id ${dir_out}/fasta ${dir_out}/bed ${dir_out}/gff
    echo "---Analysis directory (dir_out): ${dir_out}---"
    echo "------- proteome_id: ${proteome_id} -------"
    echo "------- dir_inpt: ${dir_inpt} -------"
    echo "------- flank: ${FLANK_BP} -------"
    echo "------- query: ${fa_pep_query} -------"
    echo "------- fa_pep: ${fa_pep} -------"
    echo "------- fa_assembly: ${fa_assembly} -------"
    echo "------- gff: ${gff} -------"
    echo "------- map_patterns: ${map_patterns} -------"
    ## reciprocal BLASTP
    echo "---Executing reciprocal BLASTP of query genes against proteins---"
    ${rblastgen} --dir ${dir_out}/rblastgen --seq-type p \
                 --fasta ${fa_pep_query} --bg-fasta ${fa_pep_bg} \
                 --query-fasta ${fa_pep} --relax
    ## remove blast folders because they might be massive depending on the queries
    rm -r ${dir_out}/rblastgen/blast1 ${dir_out}/rblastgen/blast2
    ## if no potential homologues, exit
    if [ ! -f ${dir_out}/rblastgen/final/recipblast.final.fasta ]; then
        return
    fi
    ## get subject ids
    echo "---Retrieving potential homologous proteins---"
    local fa_pep_rblast=${dir_out}/rblastgen/final/recipblast.final.fasta
    local txt_fid=${dir_out}/id/${proteome_id}.subset-seed.fid.txt
    local tmp=$(mktemp -p ${dir_out} tmp.XXX)
    grep '>' ${fa_pep_rblast} | sed 's/^>//' > ${tmp}
    generate_fids_v1 ${tmp} ${txt_fid} "${map_patterns}"
    rm ${tmp}
    ## subset GFF for all entries related to transcripts
    local gff_rblast=${dir_out}/gff/recipblast.gff
    ${subset_gff} -f ${txt_fid} -i ${gff} -o ${gff_rblast} --family
    ## subset genome and annotations to include some flank
    echo "---Subsetting genome and annotation around potential homologues---"
    ## get IDs of all relevant genes
    local txt_gid=${dir_out}/id/${proteome_id}.subset-seed.gid.txt
    awk '$3=="gene"' ${gff_rblast} | cut -f9 | grep -oP '^ID=[^;]+|;ID=[^;]+' | sed 's/^ID=//' > ${txt_gid}
    ## subset genome and annotations to include
    local gff_subset=${dir_out}/gff/${proteome_id}.subset.gff
    local gff_subset_og=${dir_out}/gff/${proteome_id}.subset-OGcoords.gff
    local fa_subset=${dir_out}/fasta/${proteome_id}.subset.fasta
    local bed_subset=${dir_out}/bed/${proteome_id}.subset.bed
    /mnt/chaelab/rachelle/src/extract_seq_and_annotate_merge.sh --gff ${gff} --fasta ${fa_assembly} --id ${txt_gid} \
                                                                --gff-out ${gff_subset} --fasta-out ${fa_subset} \
                                                                --bed-out ${bed_subset} \
                                                                --extend ${EXTEND_BP} --flank ${FLANK_BP} \
                                                                --gff-original ${gff_subset_og} \
                                                                --parent --memsave
    ## get genbank sequence
    if [ ${GENBANK_OUT} -eq 1 ]; then
        python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src/'); import data_manip as dm; import extract_seq_and_annotate as esa; afasta = esa.AnnotatedFasta('${fa_assembly}', '${gff_subset_og}', name = '${proteome_id}'); fids = dm.splitlines('${txt_gid}'); esa.get_flanked_seq(afasta, fids = fids, flank = ${FLANK_BP}, gb_dir = '${dir_out}/genbank', id_as_label = True)"
    fi
}

process_proteome_general_v1() {
    echo "function: process_proteome_general_v1"
    local proteome_id="${1}"
    local dir_out=${DIR}/results/${proteome_id}
    mkdir -p ${dir_out}
    ## get relevant filenames/values
    fa_assembly=$(get_genome_lookup_v1 ${proteome_id} fasta_assembly)
    fa_pep=$(get_genome_lookup_v1 ${proteome_id} fasta_protein)
    gff=$(get_genome_lookup_v1 ${proteome_id} gff)
    patterns=$(get_genome_lookup_v1 ${proteome_id} pid_gff_map)
    ## find homologues
    find_general_homologues_v1 ${proteome_id} ${dir_out} ${PEP_QUERY} ${PEP} \
                               ${fa_pep} ${fa_assembly} ${gff} "${patterns}"
    echo "=== completed proteome '${proteome_id}' ==="
}

## parse proteome IDs to process
mkdir -p ${DIR}/id
f_proteome_id=${DIR}/id/query.proteome.txt
if [ "${PROTEOME_ID}" == "all" ]; then
    grep -v '^#' ${GENOME_LOOKUP} | cut -f1 > ${f_proteome_id}
else
    echo "${PROTEOME_ID}" | tr ',' '\n' > ${f_proteome_id}
fi
proteome_ids=( $(cat ${f_proteome_id} | sort | uniq | tr '\n' ' ') )

## start processing
mkdir -p ${DIR}/results
mkdir -p ${DIR}/log
for proteome_id in ${proteome_ids[@]}; do
    (
        echo "============::=== Processing proteome '${proteome_id}' ==::============"
        logfile=${DIR}/log/${proteome_id}.log
        process_proteome_general_v1 ${proteome_id} > ${logfile} 2>&1
    ) &
    ## allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge ${PROCESSES} ]]; then
        ## now there are $N jobs already running, so wait here for any job
        ## to be finished so there is a place to start the next one
        wait -n
    fi
done
wait
echo "===::== completed ==::==="

