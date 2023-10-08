#!/bin/bash

## template paths
XML_DIR=/mnt/chaelab/rachelle/programmes/genomegaMap
XML_CONSTANT=${XML_DIR}/genomega.constant.Bayesian.xml.template
XML_CONSTANT_ML=${XML_DIR}/genomega.constant.ML.xml.template
XML_WINDOW=${XML_DIR}/genomega.window.Bayesian.xml.template
XML_SITE=${XML_DIR}/genomega.site.Bayesian.xml.template

## defaults
NITER_CONSTANT=10000
NITER_WINDOW=500000
NITER_SITE=1000000

## help message
print_help () {
    echo "Usage: run_genomega.sh [OPTIONS...] <mode> <input fasta> <output txt>"
    echo "  Output txt file not required for constantML mode"
    echo "Valid modes: constant, constantML, window/sliding, site"
    echo "OPTIONS: -i|--niter <number of iterations (for Bayesian; default: ${NITER_CONSTANT} [constant], ${NITER_WINDOW} [window/sliding], ${NITER_SITE} [site])>"
}

## parse arguments
# INPT=( "${@}" )

## shift twice after a keyword argument to skip over its value because otherwise anything that isn't a keyword argument will be added to INPT
while (( "$#" )); do
    case "$1" in
        -i|--niter) niter="${2}";shift;; ## number of iterations
        # -h|--help) man -l ${SCRIPT_DIR}/RUN_GENOMEGAMAP.1; exit 0;;
        -h|--help) print_help; exit 0;;
        *) INPT+=( "${1}" ); ## positional args
    esac
    shift
done

## check arguments
if [ ${#INPT[@]} -lt 2 ] || [ ${#INPT[@]} -gt 3 ] ; then
    echo "This script requires at least 2 and no more than 3 positional arguments"
    print_help
    exit 1
fi

## extract arguments
mode=${INPT[0]}
fasta=${INPT[1]}
output=${INPT[2]}

## select relevant template XML file & niter for user-selected mode
if [ "${mode}" == "constantML" ]; then
    xml=${XML_CONSTANT_ML}
elif [ "${mode}" == "constant" ]; then
    xml=${XML_CONSTANT}
    niter=${niter:-${NITER_CONSTANT}}
elif [ "${mode}" == "window" ] || [ "${mode}" == "sliding" ]; then
    xml=${XML_WINDOW}
    niter=${niter:-${NITER_WINDOW}}
elif [ "${mode}" == "site" ]; then
    xml=${XML_SITE}
    niter=${niter:-${NITER_SITE}}
else
    echo "Invalid mode: ${mode}"
fi

## make XML configuration file for current run
tmp=$(mktemp "genomegaMap.XXX")
cp ${xml} ${tmp}
sed -i "s|:INPUT_FILE:|${fasta}|" ${tmp}
sed -i "s|:OUTPUT_FILE:|${output}|" ${tmp}
sed -i "s|:NITER:|${niter}|" ${tmp}

## execute genomegaMap
docker run -v $(pwd):/home/ubuntu dannywilson/genomegamap ${tmp}

## remove XML configuration file
rm ${tmp}

## exit
exit 0
