#!/bin/bash

usage() {
cat << EOF
Usage: genomic_neighbourhood.sh [OPTION]...
Get overlapping and nearby features from GFF3 file for each BED range.

File control:
       -g, --gff=FILE     input GFF3 file
       -r, --range=FILE   BED file of query ranges
       -o, --out, --fout=FILE  output neighbour file (tab-separated)

Output control:
       -d, --distance=INT maximum distance (bp) in which to search for flanking
                            features (default=10000)
       -f|--flank=INT     number of features up and downstream to report
                            (default=1)
       -F|--feature=TEXT  comma-separated feature types to report
                            (default='gene,pseudogene')
       --subsetted        GFF3 file has already been subsetted and is small
                            enough to load in memory
       --sorted           GFF3 file has been sorted by 'sort -k1,1n' and there
                            is no need to re-sort it even if it exceeds the
                            chromosome length threshold
       --threshold, --sort-threshold=INT  chromosome size threshold for
                            automatic sorting of GFF3 file by 'sort -k1,1n'
                            before subsetting. Ignored if --subsetted or
                            --sorted are raised. (default=511999999)
       --dir, --directory=DIR  move intermediate files to this directory upon exit

Miscellaneous:
       -q, --quiet        silence print statments
       -v, --verbose      print information statements (default)
EOF
}


FLANK_DEFAULT=1
MAX_DISTANCE_BP_DEFAULT=10000
MEMSAVE_DEFAULT=1
SUBSETTED_DEFAULT=0
SORTED_DEFAULT=0
FEATURES_DEFAULT='gene,pseudogene'
THRESHOLD_DEFAULT=511999999 ## maximum chromosome size for bedtools without -sorted option is 512 Mb (https://bedtools.readthedocs.io/en/latest/)
QUIET_DEFAULT=0

while (( "$#" )); do
    case "$1" in
        -g|--gff) GFF="$(realpath ${2})"; shift;; ## input GFF3 file
        -r|--range) F_RANGE="$(realpath ${2})"; shift;; ## input BED file of ranges (1 per line; format: <chrom>\t<start>\t<end>)
        -o|--out|--fout) FOUT="$(realpath ${2})"; shift;; ## output file
        -d|--distance) MAX_DISTANCE_BP="${2}"; shift;; ## number of bp to extend for feature subsetting
        -f|--flank) FLANK="${2}"; shift;; ## number of features up and downstream to report
        -F|--feature) FEATURES="${2}"; shift;; ## feature types to report
        ## move intermediate files to DIR upon exit instead of deleting it
        --dir|--directory) DIR="$(realpath ${2})"; shift;;
        ## chromosome size threshold for automatically setting SORTED=1
        --sort-threshold|--threshold) THRESHOLD="${2}"; shift;;
        --sorted) SORTED=1;;
        --subsetted) SUBSETTED=1;; ## raise if GFF3 file has already been subsetted (i.e. is feasible to load in memory)
        --memsave) MEMSAVE=1;; ## use memsave mode
        --no-memsave) MEMSAVE=0;; ## don't use memsave mode
        ## execution control
        -q|--quiet) QUIET=1;; ## silence print statements
        -v|--verbose) QUIET=0;; ## print information statements
        ## print usage
        "") usage; exit 1;;
        *) usage; exit 1;;
    esac
    shift
done

## set some variables
MAX_DISTANCE_BP=${MAX_DISTANCE_BP:-${MAX_DISTANCE_BP_DEFAULT}}
FLANK=${FLANK:-${FLANK_DEFAULT}}
QUIET=${QUIET:-${QUIET_DEFAULT}}
SORTED=${SORTED:-${SORTED_DEFAULT}}
SUBSETTED=${SUBSETTED:-${SUBSETTED_DEFAULT}}
MEMSAVE=${MEMSAVE:-${MEMSAVE_DEFAULT}}
FEATURES=${FEATURES:-${FEATURES_DEFAULT}}
THRESHOLD=${THRESHOLD:-${THRESHOLD_DEFAULT}}

## define print function
function state {
    if [ ${QUIET} -eq 0 ]; then
        echo "${1}"
    fi
}

## set temporary directory
dir_tmp=$(mktemp -d)
function mktemp_f { echo "$(mktemp -p ${dir_tmp} ${1})"; }

## cleanup
function cleanup_delete { echo "$(basename $0): Removing temporary directory '${dir_tmp}'"; rm -r ${dir_tmp}; }
function cleanup_move { echo "$(basename $0): Moving temporary directory '${dir_tmp}' to '${DIR}'"; mv ${dir_tmp} ${DIR}; }
if [ -z ${DIR} ]; then
    trap cleanup_delete EXIT
    trap cleanup_delete SIGINT
else
    state "Output directory: ${DIR}"
    trap cleanup_move EXIT
    trap cleanup_move SIGINT
fi

## print some informational stuff
state "Temporary directory: '${dir_tmp}'"


#######################
##  SUBSET GFF FILE  ##
#######################

if [[ ${SUBSETTED} -ne 1 ]]; then
    ## check if GFF file needs to be sorted
    if [[ ${SORTED} -ne 1 ]]; then
        largest_chom_size=$(cut -f5 ${GFF} | sort -k1,1n | tail -1)
        if [[ ${largest_chom_size} -ge ${THRESHOLD} ]]; then
            ###########################
            echo "--Sorting GFF file--"
            ###########################
            tmp_gff_sorted=$(mktemp_f)
            grep -v '^#' ${GFF} | sort -k1,1 -k4,4n -k5,5n > ${tmp_gff_sorted}
            GFF=${tmp_gff_sorted} ## replace $GFF with a temporary file
            SORTED=1 ## set GFF file to sorted
        fi
    fi
    ###############################################
    echo "--Subsetting GFF for features in range--"
    ###############################################
    ## expand query range
    tmp_bedext=$(mktemp_f)
    awk -v f=${MAX_DISTANCE_BP} 'OFS="\t" {$2=$2-f; $3=$3+f; print $1,$2,$3}' ${F_RANGE} |
        bedtools merge -i stdin > ${tmp_bedext}
    ## intersect
    tmp_gff=$(mktemp_f)
    if [[ ${SORTED} -eq 0 ]]; then
        bedtools intersect -a ${tmp_bedext} -b ${GFF} -wb > ${tmp_gff}
    else
        bedtools intersect -a ${tmp_bedext} -b ${GFF} -wb -sorted > ${tmp_gff}
    fi
    GFF=${tmp_gff}
fi


####################
##  GET FEATURES  ##
####################

###########################
echo "--Getting features--"
###########################

# ## filter for feature types and sort the subsetted GFF file
# tmp_features=$(mktemp_f)
# tmp_gff_sorted=$(mktemp_f)
# echo "${FEATURES}" | tr ',' '\n' > ${tmp_features}
# awk -F'\t' 'NR==FNR{a[$1]; next} $3 in a' ${tmp_features} <(grep -v '^#' ${GFF}) |
#     sort -k1,1 -k4,4n -k5,5n > ${tmp_gff_sorted}

## sort the subsetted GFF file
tmp_gff_sorted=$(mktemp_f)
grep -v '^#' ${GFF} | sort -k1,1 -k4,4n -k5,5n > ${tmp_gff_sorted}

## get features
python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src/'); import genomic_neighbourhood_functions as gn; gn.get_neighbourhoods_from_files_v1(f_bed = '${F_RANGE}', f_gff = '${tmp_gff_sorted}', f_out = '${FOUT}', feature_type = '${FEATURES}'.split(','), num_flank = ${FLANK}, max_distance = ${MAX_DISTANCE_BP})"
