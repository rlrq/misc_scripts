#!/bin/bash

usage() {
cat << EOF
Usage: programme INPUT_FILE OUTPUT_FILE [OPTION]...
Convert TSV file output by BLAST (-outfmt 6) to BED file.
       INPUT_FILE ....... TSV file output by BLAST (-outfmt 6)
       OUTPUT_FILE ...... output file

Input control:
       -C, --chrom=INT    index of chromosome/molecule column (1-index)
                            (default=2 [sacc(ver)])
       --c1=INT           index of start/end position column
                            (default=9 [sstart])
       --c2=INT           index of other start/end position column
                            (default=10 [send])
       -s, --sep=STRING   separator character (default='\t')

Output control:
       -k, --keep         keep original BLAST columns after the 3 BED columns
                            (default)
       -t, --throw        discard original BLAST columns (toggle with -k/--keep)
       -S, --sort         sort output BED file by coordinates (default)
       --unsort           keep original BLAST entry order
       -d, --dir, --directory=DIR  move intermediate files to this directory upon exit

Use no more than one of -k or -t. Use no more than one of -S or --unsort.

Miscellaneous:
       -q, --quiet        silence print statments
       -v, --verbose      print information statements (default)
EOF
}

## column index for chromosome/molecule, and the 2 range values
CHROM_DEFAULT=1
COORD1_DEFAULT=9
COORD2_DEFAULT=10
SEP_DEFAULT='\t'
KEEP_DEFAULT=1
KEEP_DIR_DEFAULT=0
QUIET_DEFAULT=0
SORT_DEFAULT=1

## parse positional arguments
BLAST_IN=$1
BED_OUT=$2
if [ -z ${BLAST_IN} ]; then
    usage
    exit 1
elif [ -z ${BED_OUT} ]; then
    usage
    exit 1
fi
shift 2

## parse keyword arguments
while (( "$#" )); do
    case "$1" in
        -C|--chrom) CHROM="${2}"; shift;; ## chromosome column index
        --c1) COORD1="${2}"; shift;; ## 1st coordinate column index
        --c2) COORD2="${2}"; shift;; ## 2nd coordinate column index
        -s|--sep) SEP="${2}"; shift;; ## input file column separator character
        -k|--keep) KEEP=1;; ## keep BLAST columns
        -t|--throw) KEEP=0;; ## discard BLAST columns
        -d|--dir|--directory) DIR="${2}"; KEEP_DIR=1; shift;; ##
        -q|--quiet) QUIET=1;; ## silence print statements
        -v|--verbose) QUIET=0;; ## print information statements
        -S|--sort) SORT=1;; ## sort output BED file by coordinates
        --unsort) SORT=0;; ## keep original entry order
        "") usage; exit 1;;
        *) usage; exit 1;;
    esac
    shift
done

## set variables
CHROM=${CHROM:-${CHROM_DEFAULT}}
COORD1=${COORD1:-${COORD1_DEFAULT}}
COORD2=${COORD2:-${COORD2_DEFAULT}}
SEP=${SEP:-${SEP_DEFAULT}}
KEEP=${KEEP:-${KEEP_DEFAULT}}
KEEP_DIR=${KEEP_DIR:-${KEEP_DIR_DEFAULT}}
QUIET=${QUIET:-${QUIET_DEFAULT}}

## define print function
function state {
    if [ ${QUIET} -eq 0 ]; then
        echo "${1}"
    fi
}

## set temporary directory
dir_tmp=$(mktemp -d)
function mktemp_f { echo "$(mktemp -p ${dir_tmp})"; }

echo "${dir_tmp}"

## cleanup
function cleanup_delete { echo "Removing temporary directory '${dir_tmp}'"; rm -r ${dir_tmp}; }
echo "made cleanup_delete"
function cleanup_move { echo "Moving temporary directory '${dir_tmp}' to '${DIR}'"; mv ${dir_tmp} ${DIR}; }
echo "made cleanup_move"
if [ -z ${DIR} ]; then
    trap cleanup_delete EXIT
    trap cleanup_delete SIGINT
else
    trap cleanup_move EXIT
    trap cleanup_move SIGINT
fi

## print some informational stuff
state "Intermediate files will be written to '${DIR}'"

## process
unsorted_bed=${dir_tmp}/unsorted.bed
if [ ${KEEP} -eq 1 ]; then
    awk -v sep=${SEP} -v C=${CHROM} -v c1=${COORD1} -v c2=${COORD2} 'BEGIN {FS=sep; OFS=sep} {min=$c1;max=$c2;if($c2<$c1)min=$c2;if($c2<$c1)max=$c1;print $C,min-1,max,$0}' ${BLAST_IN} > ${unsorted_bed}
else
    awk -v sep=${SEP} -v C=${CHROM} -v c1=${COORD1} -v c2=${COORD2} 'BEGIN {FS=sep; OFS=sep} {min=$c1;max=$c2;if($c2<$c1)min=$c2;if($c2<$c1)max=$c1;print $C,min-1,max}' ${BLAST_IN} > ${unsorted_bed}
fi

## sort
if [ ${SORT} -eq 1 ]; then
    sort -k1,1 -k2,2n, -k3,3n ${unsorted_bed} > ${BED_OUT}
else
    cp ${unsorted_bed} ${BED_OUT}
fi

exit 0
