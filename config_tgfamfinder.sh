#!/bin/bash

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
DIR_BIN_TGFAM=/mnt/chaelab/rachelle/programmes/TGFam-Finder/scripts_ubuntu/scripts/

## DEFAULTS
DEFAULT_PROGRAM_PATH_TEMPLATE=/mnt/chaelab/rachelle/programmes/TGFam-Finder/PROGRAM_PATH.config.template
DEFAULT_RESOURCE_TEMPLATE=/mnt/chaelab/rachelle/programmes/TGFam-Finder/RESOURCE.config.template
DEFAULT_EXTENSION_LENGTH='10000'
DEFAULT_MAX_INTRON_LENGTH='10000'
DEFAULT_THREADS='1'

## DEFAULTS: OPTIONAL INFORMATION
# DEFAULT_CDS_OF_TARGET_GENOME=''
# DEFAULT_GFF3_OF_TARGET_GENOME=''
# DEFAULT_RNASEQ_REVERSE_PATH=''
# DEFAULT_EXCLUDED_DOMAIN_ID=''
# DEFAULT_HMM_MATRIX_NAME=''
DEFAULT_HMM_CUTOFF='1e-3'

# params="$@"

while (( "$#" )); do
    case "$1" in
        --program-path-template) PROGRAM_PATH_TEMPLATE="$(realpath "${2}")";;
        --resource-template) RESOURCE_TEMPLATE="$(realpath "${2}")";;
        --running-path) RUNNING_PATH="$(realpath "${2}")";;
        --target-genome) TARGET_GENOME="$(realpath "${2}")";;
        --proteins-for-domain-identification) PROTEINS_FOR_DOMAIN_IDENTIFICATION="$(realpath "${2}")";;
        --tsv-for-domain-identification) TSV_FOR_DOMAIN_IDENTIFICATION="$(realpath "${2}")";;
        --resource-protein) RESOURCE_PROTEIN="$(realpath "${2}")";;
        --blast-db-name) BLAST_DB_NAME="$(realpath "${2}")";;
        --output-prefix) OUTPUT_PREFIX="${2}";;
        --target-domain-id) TARGET_DOMAIN_ID="${2}";;
        --target-domain-name) TARGET_DOMAIN_NAME="${2}";;
        --representative-domain-name) REPRESENTATIVE_DOMAIN_NAME="${2}";;
        --extension-length) EXTENSION_LENGTH="${2}";;
        --max-intron-length) MAX_INTRON_LENGTH="${2}";;
        --threads) THREADS="${2}";;
        --cds-of-target-genome) CDS_OF_TARGET_GENOME="$(realpath "${2}")";;
        --gff3-of-target-genome) GFF3_OF_TARGET_GENOME="$(realpath "${2}")";;
        --rnaseq-forward-path) RNASEQ_FORWARD_PATH="$(realpath "${2}")";;
        --rnaseq-reverse-path) RNASEQ_REVERSE_PATH="$(realpath "${2}")";;
        --excluded-domain-id) EXCLUDED_DOMAIN_ID="${2}";;
        --hmm-matrix-name) HMM_MATRIX_NAME="$(realpath "${2}")";;
        --hmm-cutoff) HMM_CUTOFF="${2}";;
        --overwrite) OVERWRITE="T";;
    esac
    shift
done

## check that required args are provided
missing_msg() {
    arg=${1}
    descr=${2}
    echo "Missing argument: ${arg} <${descr}>"
}

if [ -z "${RUNNING_PATH}" ]; then
    missing_msg "--running-path" "Location of directory running TGFam-Finder"
    exit 1
elif [ -z "${TARGET_GENOME}" ]; then
    missing_msg "--target-genome" "Full path of assembled genome (ex, path/filename), fasta format"
    exit 1
elif [ -z "${PROTEINS_FOR_DOMAIN_IDENTIFICATION}" ]; then
    missing_msg "--proteins-for-domain-identification" "Full path of peptide sequences of target or allied species, fasta format"
    exit 1
elif [ -z "${TSV_FOR_DOMAIN_IDENTIFICATION}" ]; then
    missing_msg "--tsv-for-domain-identification" "Full path of InterPro result of the peptide sequences, tsv format"
    exit 1
elif [ -z "${RESOURCE_PROTEIN}" ]; then
    missing_msg "--resource-protein" "Full path of target peptide sequences in multiple species that users prepare as resources for each annotation step, fasta format"
    exit 1
elif [ -z "${OUTPUT_PREFIX}" ]; then
    missing_msg "--output-prefix" "Output prefix generated results"
    exit 1
elif [ -z "${TARGET_DOMAIN_ID}" ]; then
    missing_msg "--target-domain-id" "Target domain ID saved in fifth column of tsv files such as PF00000, SSF00000, and SM00000, not IPR number; comma-separated if multiple (e.g. 'PF00000,SSF00000')"
    exit 1
elif [ -z "${TARGET_DOMAIN_NAME}" ]; then
    missing_msg "--target-domain-name" "target domain name(s); comma-separated if multiple (e.g. NB-ARC,TIR), same order as --target-domain-id"
    exit 1
elif [ -z "${REPRESENTATIVE_DOMAIN_NAME}" ]; then
    missing_msg "--representative-domain-name" "Target gene-family name"
    exit 1
fi

## set default values for optional args
if [ -z "${BLAST_DB_NAME}" ]; then
    missing_msg "--blast-db-name" "Full path of BLAST DB that will be automatically created"
    BLAST_DB_NAME="${RUNNING_PATH}/db/db"
    printf "\tSetting default DB location: ${BLAST_DB_NAME}\n"
fi
PROGRAM_PATH_TEMPLATE="${PROGRAM_PATH_TEMPLATE:-${DEFAULT_PROGRAM_PATH_TEMPLATE}}"
RESOURCE_TEMPLATE="${RESOURCE_TEMPLATE:-${DEFAULT_RESOURCE_TEMPLATE}}"
EXTENSION_LENGTH="${EXTENSION_LENGTH:-${DEFAULT_EXTENSION_LENGTH}}"
MAX_INTRON_LENGTH="${MAX_INTRON_LENGTH:-${DEFAULT_MAX_INTRON_LENGTH}}"
THREADS="${THREADS:-${DEFAULT_THREADS}}"
HMM_CUTOFF="${HMM_CUTOFF:-${DEFAULT_HMM_CUTOFF}}"

## make .config file paths in ${RUNNING_PATH} and copy templates to them
if [ -d "${RUNNING_PATH}" ] && ! [ -z "${OVERWRITE}" ]; then
    echo "${RUNNING_PATH} already exists. Programme aborted. To allow overwriting of existing directories, use '--overwrite'."
    exit 2
fi
mkdir -p "${RUNNING_PATH}"
RESOURCE_CONFIG="${RUNNING_PATH}/RESOURCE.config"
PROGRAM_PATH_CONFIG="${RUNNING_PATH}/PROGRAM_PATH.config"
cp "${RESOURCE_TEMPLATE}" "${RESOURCE_CONFIG}"
cp "${PROGRAM_PATH_TEMPLATE}" "${PROGRAM_PATH_CONFIG}"

## set params
PROGRAM_PATH_PARAMS=( RUNNING_PATH )
RESOURCE_PARAMS=( TARGET_GENOME PROTEINS_FOR_DOMAIN_IDENTIFICATION TSV_FOR_DOMAIN_IDENTIFICATION RESOURCE_PROTEIN BLAST_DB_NAME OUTPUT_PREFIX TARGET_DOMAIN_ID TARGET_DOMAIN_NAME REPRESENTATIVE_DOMAIN_NAME EXTENSION_LENGTH MAX_INTRON_LENGTH THREADS CDS_OF_TARGET_GENOME GFF3_OF_TARGET_GENOME RNASEQ_FORWARD_PATH RNASEQ_REVERSE_PATH EXCLUDED_DOMAIN_ID HMM_MATRIX_NAME HMM_CUTOFF )

for param in ${PROGRAM_PATH_PARAMS[@]}; do
    sed -i "s|:${param}:|${!param}|" "${PROGRAM_PATH_CONFIG}"
done

for param in ${RESOURCE_PARAMS[@]}; do
    sed -i "s|:${param}:|${!param}|" "${RESOURCE_CONFIG}"
done

## execute Build_TGFam-Finder.pl
cd "${RUNNING_PATH}"
perl "${DIR_BIN_TGFAM}/Build_TGFam-Finder_rlrq20230623.pl" "${RUNNING_PATH}" RESOURCE.config PROGRAM_PATH.config
