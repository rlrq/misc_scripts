#!/bin/bash

## TODO: accommodate GFF3 instead of requiring the GFF2BED version
ORIGINAL_DIR=$(pwd)
# SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

REF_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta'
BED_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed'
GFF_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.gff'
TOPN_DEFAULT=1
METRIC_DEFAULT='bitscore'
GETSEQ_DEFAULT='/mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh'

while (( "$#" )); do
    case "$1" in
        ## input options
        -q|--query) QUERY="${2}";;
        -r|--ref|--reference) REF="${2}";; ## fasta file
        -b|--bed) BED="${2}";;
        -g|--gff) GFF="${2}";;
        ## execution options
        -n|--topn) TOPN="${2}";;
        -m|--metric) METRIC="${2}";;
        ## output options
        --fout-blast) FOUT_BLAST="${2}";;
        --fout-bed) FOUT_BED="${2}";;
        --fout-gff) FOUT_GFF="${2}";;
        --fout-gid) FOUT_GID="${2}";;
        --fout-intersect) FOUT_INTERSECT="${2}";;
        # --fout-intersect-full) FOUT_FULL="${2}";;
        --fout-gene) FOUT_GENE="${2}";;
        --fout-cds) FOUT_CDS="${2}";;
        --fout-protein) FOUT_PEP="${2}";;
        --fout-complete) FOUT_COMPLETE="${2}";;
        --getseq) GETSEQ="${2}";;
        # --fmt) FMT="${2}";;
        # --acc-lookup) ACC_LOOKUP="${2}"; MODE='ref';;
        # -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        # --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

DIR="${DIR:-${DIR_DEFAULT}}"
REF="${REF:-${REF_DEFAULT}}"
BED="${BED:-${BED_DEFAULT}}"
GFF="${GFF:-${GFF_DEFAULT}}"
TOPN="${TOPN:-${TOPN_DEFAULT}}"
METRIC="${METRIC:-${METRIC_DEFAULT}}"
GETSEQ="${GETSEQ:-${GETSEQ_DEFAULT}}"

## temporary directory stuff
tmpdir=$(mktemp -d)
echo "tmpdir: ${tmpdir}"
trap 'rm -r "${tmpdir}"' EXIT
trap 'rm -r "${tmpdir}"' SIGINT

## get gene-only annotations (format: <chrom>\t<start>\t<end>\t<gid>)
tmp_gene_bed=$(mktemp -p ${tmpdir})
if [ "${GFF}" != "${GFF_DEFAULT}" ]; then
    awk '$3=="gene" || $3=="pseudogene"' ${GFF} | gff2bed | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}'> ${tmp_gene_bed}
else
    awk 'BEGIN{OFS="\t"} $8=="gene" || $8=="pseudogene" {print $1,$2,$3,$4}' ${BED} > ${tmp_gene_bed}
fi

## BLAST and get top n genes
source /mnt/chaelab/rachelle/src/run_blast6.sh
topn_gid () {
    local query=${1}
    local fout_intersect=${2}
    local fout_gid=${3}
    local tmp_blast=${FOUT_BLAST:-$(mktemp -p ${tmpdir})}
    local tmp_blast_bed=$(mktemp -p ${tmpdir})
    ## blast to reference
    echo "--BLAST--"
    run_blast6 blastn -query ${query} -subject ${REF} -out ${tmp_blast} \
               -outfmt "6 qseqid qstart qend sseqid sstart send sstrand ${METRIC}"
    ## get top n hits, coerce into BED
    ## ## didn't move this section to R script file in order to keep script to 1 file; side effect: it's printed
    echo "--Getting top N hits--"
    R -q -e "suppressMessages(suppressWarnings(library(tidyverse))); df <- read.table('${tmp_blast}', header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>% group_by(qseqid) %>% top_n(${METRIC}, n = ${TOPN}) %>% ungroup() %>% dplyr::select(sseqid, sstart, send, qseqid, qstart, qend, ${METRIC}); df <- df %>% dplyr::mutate(start = sapply(1:nrow(df), function(i){min(df[i, c('sstart', 'send')]) - 1}), end = sapply(1:nrow(df), function(i){max(df[i, c('sstart', 'send')])})) %>% dplyr::select(sseqid, start, end, qseqid, qstart, qend, ${METRIC}) %>% arrange() %>% write.table('${tmp_blast_bed}', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')"
    ## intersect w/ BED/GFF
    echo "--Intersecting w/ reference annotations--"
    if [ -z "${BED_ARGS}" ]; then
        bedtools_options="-wao"
    else
        bedtools_options="-wao ${BED_ARGS}"
    fi
    bedtools intersect ${bedtools_options} -a ${tmp_blast_bed} -b ${tmp_gene_bed} > ${fout_intersect}
    ## write GID
    cut -f11 ${fout_intersect} | sort | uniq | grep -v '\.' > ${fout_gid}
    ## remove tmp files
    rm ${tmp_blast_bed}
}

## execute
tmp_gid=$(mktemp -p ${tmpdir})
if (! [ -z "${FOUT_INTERSECT}" ]); then
    topn_gid ${QUERY} ${FOUT_INTERSECT} ${tmp_gid}
else
    topn_gid ${QUERY} $(mktemp -p ${tmpdir}) ${tmp_gid}
fi

## write FOUT_GID if requested
if (! [ -z "${FOUT_GID}" ]); then
    cp ${tmp_gid} ${FOUT_GID}
fi

## write FOUT (GFF/BED)
tmp_fout_gff=${FOUT_GFF:-$(mktemp -p ${tmpdir})}
tmp_fout_bed=${FOUT_BED:-$(mktemp -p ${tmpdir})}
if (! [ -z "${FOUT_GFF}" ] || ! [ -z "${FOUT_BED}" ] || \
   ! [ -z "${FOUT_GENE}" ] || ! [ -z "${FOUT_CDS}" ] || ! [ -z "${FOUT_PEP}" ]); then
    echo "--Writing annotations--"
    ## make gff file
    if [ "${GFF}" != "${GFF_DEFAULT}" ]; then
        /mnt/chaelab/rachelle/src/extract_gff_features.py ${GFF} ${tmp_gid} ${tmp_fout_gff}
    else
        /mnt/chaelab/rachelle/src/extract_gff_features.py ${BED} ${tmp_gid} ${tmp_fout_gff}
    fi
    ## make bed file
    /mnt/chaelab/rachelle/src/gff_to_bed.sh ${tmp_fout_gff} ${tmp_fout_bed}
fi

## write sequences
if (! [ -z "${FOUT_GENE}" ] || ! [ -z "${FOUT_CDS}" ] || ! [ -z "${FOUT_PEP}" ]); then
    echo "--Writing sequences--"
    getseq_common="${GETSEQ} -g ${tmp_gid} -a ref -r ${REF} -b ${tmp_fout_bed} --adjust-dir --no-bed"
    if (! [ -z "${FOUT_GENE}" ]); then
        ${getseq_common} -f gene -o ${FOUT_GENE}
    fi
    if (! [ -z "${FOUT_CDS}" ]); then
        ${getseq_common} -f CDS -o ${FOUT_CDS}
    fi
    if (! [ -z "${FOUT_PROTEIN}" ]); then
        ${getseq_common} -f CDS -o ${FOUT_PEP} --translate
    fi
    if (! [ -z "${FOUT_COMPLETE}" ]); then
        ${getseq_common} -f CDS -o ${FOUT_COMPLETE} --complete
    fi
fi

