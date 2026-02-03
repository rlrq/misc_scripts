#!/bin/bash

usage() {
cat << EOF
Usage: getcol.sh [OPTION] FILE
Subset a delimited file by column name.

FILE   input delimited file

File control:
       -F, --file=FILE    file of column names (one per line)
       -c|--col|--cols|--colname|--colnames=STRING  comma-separated column names
       -d|--delim|-FS|--FS=STRING  delimiter of input file (default='\t')

-F and -c are mutually exclusive. Output will be printed to stdout.

EOF
}

FS_DEFAULT='\t'
OFS_DEFAULT='\t'

while (( "$#" )); do
    case "$1" in
        -F|--file) FILE="${2}"; shift;; ## file of column names (mutually exclusive with --col)
        -c|--col|--cols|--colname|--colnames) COLNAMES="${2}"; shift;; ## comma-separated colnames (mutually exclusive with -F)
        -d|--delim|-FS|--FS) FS="${2}"; shift;; ## delimiter of input file
        # --od|-OFS|--OFS) OFS="${2}"; shift;; ## delimiter of output
        *) if ! [ -z ${INPT} ]; then usage; exit 1; else INPT="${1}"; fi;; ## input file
        "") usage; exit 1;;
    esac
    shift
done

## set some variables
FS="${FS:-${FS_DEFAULT}}"
OFS="${OFS:-${OFS_DEFAULT}}"

## convert input colnames to file
f_tmp=$(mktemp)
if [ -z ${FILE} ]; then
    echo "${COLNAMES}" | tr ',' '\n' > ${f_tmp}
    FILE=${f_tmp}
fi

## subset by column and print to stdout
cut -d "${FS}" -f"$(head -n1 ${INPT} | sed "s/${FS}/\t/g" | tr '\t' '\n' | grep -Fxn -f "${FILE}" | cut -f1 -d: | paste -sd,)" "${INPT}"

## remove temporary file
rm ${f_tmp}

