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
DOMAIN_NAME_DEFAULT=domain
PREFIX_DEFAULT='get_domain'

params="$@"

while (( "$#" )); do
    case "$1" in
        -f|--fasta|-i|--input) FASTA="${2}";;
        -d|--dir) DIR="${2}";;
        --domain) DOMAIN="${2}";;
        --tsv-domain) TSV_DOMAIN="${2}";;
        --domain-name) DOMAIN_NAME="${2}";;
        --domain-db) RPS_DB="${2}";;
        --prefix) PREFIX="${2}";;
        --keep-tmp) KEEP_TMP='True';;
        -o|--out) FOUT="$(realpath ${2})";;
    esac
    shift
done

DIR="${DIR:-${DIR_DEFAULT}}"
RPS_DB="${RPS_DB:-${RPS_DB_DEFAULT}}"
KEEP_TMP="${KEEP_TMP:-False}"
DOMAIN_NAME="${DOMAIN_NAME:-${DOMAIN_NAME_DEFAULT}}"
PREFIX="${PREFIX:-${PREFIX_DEFAULT}}"

mkdir -p ${DIR}

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
    DOMAIN_NAME=${DOMAIN}
fi

## parse genes
if ! [ -z "${INPUT_DIR}" ]; then
    parsed_genes=( $(ls ${INPUT_DIR} | sed 's/\.[^\.]*//g') )
elif [ -f "${GENE}" ]; then ## if in file
    readarray -d ',' -t parsed_genes <<< "$(tr '\n' ',' < ${GENE} | sed 's/,\+$//g')"
else ## if comma-separated
    readarray -d , -t parsed_genes <<< ${GENE}
fi

if [ -z "${TSV_DOMAIN}" ]; then
    ## identify domain positions
    echo "Getting domain ranges"
    tsv_domain=${DIR}/${PREFIX}.${DOMAIN}.tsv
    tmp_f=${DIR}/${PREFIX}.rpsblast.tsv
    tmp_f2=${DIR}/${PREFIX}.rpsblast.tsv.tmp
    /mnt/chaelab/rachelle/programmes/BLAST+/ncbi-blast-2.11.0+/bin/rpsblast -db ${RPS_DB} -query ${FASTA} -out ${tmp_f} \
                                                                            -outfmt "6 $(echo ${HEADER_CDD} | tr ',' ' ')"
    ## generate data for input as "DOMAIN_F" to getSeq
    awk -v raw="${domain}" -v d="${DOMAIN}" '{if ($2 ~ raw) print $0 "\t" d}' ${tmp_f} >> ${tmp_f2}
    ## combine tmp_f into final DOMAIN_F for getSeq
    echo -e "$(echo ${HEADER_CDD},domain | tr ',' '\t')\n$(cat ${tmp_f2})" > ${tsv_domain}
    TSV_DOMAIN=${tsv_domain}
    ## add header to tmp_f (tmp_f won't be used anymore in this run but it can be reused by future runs)
    echo -e "$(echo ${HEADER_CDD},domain | tr ',' '\t')\n$(cat ${tmp_f})" > ${tmp_f2}
    mv ${tmp_f2} ${tmp_f}
else
    echo "Getting domain ranges"
    tsv_domain=${DIR}/${PREFIX}.${DOMAIN}.tsv
    tmp_f2=${DIR}/${PREFIX}.rpsblast.tsv.tmp
    ## generate data for input as "DOMAIN_F" to getSeq
    awk -v raw="${domain}" -v d="${DOMAIN}" '{if ($2 ~ raw) print $0 "\t" d}' ${TSV_DOMAIN} >> ${tmp_f2}
    ## combine tmp_f into final DOMAIN_F for getSeq
    echo -e "$(head -1 ${TSV_DOMAIN})$(cat ${tmp_f2})" > ${tsv_domain}
    TSV_DOMAIN=${tsv_domain}
fi

## get desired sequence
echo "Getting domain sequences"
FOUT="${FOUT:-${DIR}/${PREFIX}.${DOMAIN_NAME}.fasta}"
## proceed only if >=1 entry in ${tsv_domain}
if [ $(wc -l ${TSV_DOMAIN} | grep -o "^[0-9]\+") -gt 1 ]; then
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); import data_manip, basics, fasta_manip; from Bio import SeqIO; get, data = data_manip.parse_get_data('${TSV_DOMAIN}'); domain = '${domain}'; data = [x for x in data if get(x, 'sseqid').split('|')[-1] == domain]; by_seq = {seqid: [tuple(map(int, get(x, 'qstart', 'qend'))) for x in data if get(x, 'qseqid') == seqid] for seqid in set(get(data, 'qseqid'))}; by_seq = {seqid: sorted(basics.merge_overlapping_ranges(*[(r[0]-1, r[1]) for r in ranges], end_inclusive = False)) for seqid, ranges in by_seq.items()}; records = SeqIO.parse('${FASTA}', 'fasta'); to_write = {(record.id + '|${DOMAIN_NAME}|' + str(i+1)): record.seq[r[0]:r[1]] for record in records for i, r in enumerate(by_seq.get(record.id, []))}; fasta_manip.dict_to_fasta(to_write, '${FOUT}')"
else
    echo "Domain not found"
fi

if [[ ! "${KEEP_TMP}" == "True" ]]; then
    if ! [ -z ${tmp_f} ]; then
        rm ${tmp_f}
    fi
    if [ -f ${tmp_f2} ]; then
        rm ${tmp_f2}
    fi
    if ! [ -z ${tsv_domain} ]; then
        rm ${tsv_domain}
    fi
else
    if [ -f ${tmp_f2} ]; then
        rm ${tmp_f2}
    fi
fi

exit 0
