#!/bin/bash

## NOTE: This script does not support --input-dir

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

## because aliases don't really work with non-interactive shell :(
get_seq=/mnt/chaelab/rachelle/scripts/get_seqs/get_seqs.sh

declare -A DOMAINS
DOMAINS["TIR"]='307630';DOMAINS["RX-CC_like"]='271353';DOMAINS["CC"]='271353';DOMAINS["RPW8"]='310336';DOMAINS["NB-ARC"]='307194';DOMAINS["NBS"]='307194'

HEADER_CDD="qseqid,sseqid,pident,length,qstart,qend"
DB_CDD='/mnt/chaelab/shared/blastdb/RPSdb/Cdd'
RPS_DB_DEFAULT=${DB_CDD}
DIR_DEFAULT="/mnt/chaelab/$(whoami)/get_domain"

params="$@"

while (( "$#" )); do
    case "$1" in
        -g|--gene) GENE="${2}";;
        -d|--dir) DIR="${2}";;
        --domain) DOMAIN="${2}";;
        --domain-db) RPS_DB="${2}";;
        --input-dir) INPUT_DIR="${2}";;
        --prefix) PREFIX="${2}";;
        --keep-tmp) KEEP_TMP='True';;
    esac
    shift
done

DIR="${DIR:-${DIR_DEFAULT}}"
RPS_DB="${RPS_DB:-${RPS_DB_DEFAULT}}"
KEEP_TMP="${KEEP_TMP:-False}"

mkdir -p ${DIR}
tmp_f=${DIR}/tmp.txt

if [[ ! " ${!DOMAINS[@]} " =~ " ${DOMAIN} " ]]; then
    echo "${DOMAIN} is not a supported domain name."
    if [[ ${DOMAIN} =~ ^[0-9]+$ ]]; then
        echo "Attempting to parse as CDD PSSM-Id."
        domain=${DOMAIN}
    else
        echo "Unexpected domain input. Please provide a CDD PSSM-Id (numeric) or one of the following supported domain names: $(echo ${!DOMAINS[@]} | sed 's/ /, /g')"
        exit 1
    fi
else
    domain=${DOMAINS[${DOMAIN}]}
fi


## parse genes
if ! [ -z "${INPUT_DIR}" ]; then
    parsed_genes=( $(ls ${INPUT_DIR} | sed 's/\.[^\.]*//g') )
elif [ -f "${GENE}" ]; then ## if in file
    readarray -d ',' -t parsed_genes <<< "$(tr '\n' ',' < ${GENE} | sed 's/,\+$//g')"
else ## if comma-separated
    readarray -d , -t parsed_genes <<< ${GENE}
fi

genes=$(echo ${parsed_genes[@]} | sed 's/ /,/g' | sed 's/\n//g')
echo "genes: ${genes[@]}"

## get protein sequences
echo "Getting reference protein sequences"
tsv_domain=${DIR}/${PREFIX}_${DOMAIN}.tsv
fa_aa=${DIR}/${PREFIX}_${DOMAIN}_protein.fasta
$get_seq -g ${genes} -a ref -f CDS -d ${DIR} -o ${fa_aa} --translate --no-bed > /dev/null

## identify domain positions
echo "Getting domain ranges"
rpsblast+ -db ${RPS_DB} -query ${fa_aa} -out ${tsv_domain} \
          -outfmt "6 $(echo ${HEADER_CDD} | tr ',' ' ')"
## generate data for input as "DOMAIN_F" to getSeq
awk -v raw="${domain}" -v d="${DOMAIN}" '{if ($2 ~ raw) print $0 "\t" d}' ${tsv_domain} >> ${tmp_f}
## combine tmp_f into final DOMAIN_F for getSeq
echo -e "$(echo ${HEADER_CDD},domain | tr ',' '\t')\n$(cat ${tmp_f})" > ${tsv_domain}

## get desired sequence
echo "Getting domain sequences"
## proceed only if >=1 entry in ${tsv_domain}
if [ $(wc -l ${tsv_domain} | grep -o "^[0-9]\+") -gt 1 ]; then
    $get_seq $params --domain-file ${tsv_domain} --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" > /dev/null
else
    echo "Domain not found"
fi

if [[ ! "${KEEP_TMP}" == "True" ]]; then
    rm ${tmp_f} ${tsv_domain} ${fa_aa}
else
    rm ${tmp_f}
fi

exit 0
