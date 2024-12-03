#!/bin/bash

while (( "$#" )); do
    case "$1" in
        -d|--dir) DIR="${2}"; shift;; ## output directory
        --pref|--prefix) PREFIX="${2}"; shift;; ## output prefix, used if --aln-out and/or --tree-out are not specified
        --tsv) TSV_DOMAINS="${2}"; shift;; ## one of the outputs of rpsblast_and_get_nlr_architecture.sh
        --pep|--fasta) FA_PEP="${2}"; shift;; ## fasta file of peptides to get domains for
        --domain|--domains) DOMAINS="${2}"; shift;; ## comma-separated str that will be parsed into an array
        --fasta-cds) FA_CDS="${2}"; shift;; ## optional
    esac
    shift
done

if [[ -z ${TSV_DOMAINS} ]]; then
    echo "TSV file of modified RPS-blast output required. Please provide using '--tsv <file>'. See rpsblast_and_get_nlr_architecture.sh's ...domain.tsv output file for format."
    exit 1
elif [[ -z ${FA_PEP} ]]; then
    echo "Fasta file of peptide sequences required. Please provide using '--fasta <file>' or '--pep <file>'."
    exit 1
elif [[ -z ${DOMAINS} ]]; then
    echo "Domain(s) required. Please provide using '--domain <comma-separated domain short names>'."
    exit 1
fi

DIR=${DIR:-$(pwd)}
PREFIX=${PREFIX:-extractDomains}

extract_domain() {
    local outdir=${1}
    local prefix=${2}
    local tsv_domains=${3}
    local fasta=${4}
    local domain=${5}
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from data_manip import parse_get_data; from fasta_manip import fasta_to_dict, dict_to_fasta; from range_manip import ranges_union; get, dat = parse_get_data('${tsv_domains}'); dat = [x for x in dat if get(x, 'domain')=='${domain}']; ranges = {qseqid: sorted(ranges_union([[(get(x, 'qstart') - 1, get(x, 'qend'))] for x in dat if get(x, 'qseqid') == qseqid])) for qseqid in set(get(dat, 'qseqid'))}; seqs = fasta_to_dict('${fasta}'); to_write = {qseqid + '|' + '${domain}' + '|' + str(i+1) + '|' + str(r[0]+1) + '-' + str(r[1]): seqs[qseqid][r[0]:r[1]] for qseqid, qseq_dat in ranges.items() for i, r in enumerate(qseq_dat)}; dict_to_fasta(to_write, '${outdir}/${prefix}.${domain}.fasta')"
}

extract_domain_cds () {
    local outdir=${1}
    local prefix=${2}
    local tsv_domains=${3}
    local fasta=${4} ## CDS
    local domain=${5}
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from data_manip import parse_get_data; from fasta_manip import fasta_to_dict, dict_to_fasta; from range_manip import ranges_union; get, dat = parse_get_data('${tsv_domains}'); dat = [x for x in dat if get(x, 'domain')=='${domain}']; ranges = {qseqid: sorted(ranges_union([[(get(x, 'qstart') - 1, get(x, 'qend'))] for x in dat if get(x, 'qseqid') == qseqid])) for qseqid in set(get(dat, 'qseqid'))}; seqs = fasta_to_dict('${fasta}'); to_write = {qseqid + '|' + '${domain}' + '|' + str(i+1) + '|' + str((r[0]*3)+1) + '-' + str(r[1]*3): seqs[qseqid][(r[0]*3):(r[1]*3)] for qseqid, qseq_dat in ranges.items() for i, r in enumerate(qseq_dat)}; dict_to_fasta(to_write, '${outdir}/${prefix}.${domain}.CDS.fasta')"
}

extract_domains () {
    local outdir=${1}
    local prefix=${2}
    local tsv_domains=${3}
    local fasta=${4}
    local domains=${5} ## space-delimited str that will be parsed into an array
    local fasta_cds=${6}
    for domain in ${domains[@]}; do
        extract_domain ${outdir} ${prefix} ${tsv_domains} ${fasta} ${domain}
        if ! [ -z ${fasta_cds} ]; then
            mv ${outdir}/${prefix}.${domain}.fasta ${outdir}/${prefix}.${domain}.pep.fasta
            extract_domain_cds ${outdir} ${prefix} ${tsv_domains} ${fasta_cds} ${domain}
        fi
    done
}

extract_domains ${DIR} ${PREFIX} ${TSV_DOMAINS} ${FA_PEP} "$(printf "${DOMAINS}" | tr ',' ' ')" ${FA_CDS}
