#!/bin/bash

usage() {
cat << EOF
Usage: subset_gff.sh [OPTION]...
Subset a GFF3 file by feature ID.

File control:
       -i, --in=FILE      input GFF3 file
       -o, --out=FILE     output GFF3 file
       -f, --fid=FILE     TXT file of query feature IDs (one per line)

Output control:
       -s, --self         output query features (default)
       -p, --parent       output direct parent(s) of query features
       -P, --parents, -A, --ancestors  output direct AND indirect
                            parent(s) of query features
       -c, --child        output child(ren) of query features
       -C, --children, -D, --descendants  output direct AND indirect
                            child(ren) of query features
       -F, --family       output query features, ancestors of query features AND
                            the descendants of those ancestors, regardless of if
                            they are descended from query features themselves
       --sort             sort output by chromosome, start, and end positions in
                            ascending order (sort -k1,1 -k4,4n -k5,5n)
       --unsort           keep original entry order (default)
       -d, --dir, --directory=DIR  move intermediate files to this directory upon exit

Use no more than one of -p or -P. Use no more than one of -c or -C.
Use of -F should not be accompanied by any of -s, -p, -P, -c, or -C.
Use no more than one of --sort or --unsort.

Miscellaneous:
       -q, --quiet        silence print statments
       -v, --verbose      print information statements (default)
EOF
}

## column index for chromosome/molecule, and the 2 range values
QUIET_DEFAULT=0
SORT_DEFAULT=0
CHUNK_DEFAULT=300

## parse keyword arguments
while (( "$#" )); do
    case "$1" in
        -i|--in) GFF_IN="$(realpath ${2})"; shift;; ## input GFF3 file
        -o|--out) GFF_OUT="$(realpath ${2})"; shift;; ## output GFF3 file
        -f|--fid) FID="$(realpath ${2})"; shift;; ## input feature ID TXT file (one ID per line)
        -s|--self) SELF=1;; ## output entries for query features
        -p|--parent) PARENT=1;; ## output entries for parents of query features
        -P|--parents|-A|--ancestors) ANCESTORS=1;; ## output all ancestors of query features
        -c|--child) CHILD=1;; ## output entries for children of query features
        -C|--children|-D|--descendants) DESCENDANTS=1;; ## ouptut all descendants of query features
        -F|--family) FAMILY=1;; ## output entries of ancestors + descendants + self of query features
        -d|--dir|--directory) DIR="$(realpath ${2})"; KEEP_DIR=1; shift;; ##
        -q|--quiet) QUIET=1;; ## silence print statements
        -v|--verbose) QUIET=0;; ## print information statements
        --sort) SORT=1;; ## sort output using 'sort -k1,1 -k4,4n -k5,5n'
        --unsort) SORT=0;; ## retain original entry order
        --chunk) CHUNK="${2}"; shift;; ## limit the number of feature IDs to grep at once
        "") usage; exit 1;;
        *) usage; exit 1;;
    esac
    shift
done

## set some variables
QUIET=${QUIET:-${QUIET_DEFAULT}}
SORT=${SORT:-${SORT_DEFAULT}}
CHUNK=${CHUNK:-${CHUNK_DEFAULT}}

## raise flag for descendants, ancestors, and self if --family
PARENT=${PARENT:-0}
CHILD=${CHILD:-0}
FAMILY=${FAMILY:-0}
DESCENDANTS=${DESCENDANTS:-${FAMILY}}
ANCESTORS=${ANCESTORS:-${FAMILY}}
SELF=${SELF:-${FAMILY}}

## use default mode (--self) if none of the other options are specified
if [ -z ${ANCESTORS} ] && [ -z ${DESCENDANTS} ] && [ -z ${PARENT} ] && [ -z ${CHILD} ]; then
    SELF=1
fi

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

function subset_self() {
    state "---subset_self---"
    local gff fid fout dir_self
    dir_self=${dir_tmp}
    while (( "$#" )); do
        case "$1" in
            -g|--gff) gff="${2}"; shift;;
            -f|--fid) fid="${2}"; shift;;
            -o|--out|--fout) fout="${2}"; shift;;
            -d|--dirself) dir_self="${2}"; shift;;
        esac
        shift
    done
    mkdir -p ${dir_self}
    local f_pattern=$(mktemp -p ${dir_self} self.pattern.XXXXXX)
    ## TODO: grep -fF file with ID=<fid> first, then do the Perl thing.
    sort ${fid} | uniq | ## remove duplicates
        perl -pe 's/(?=[\@#\$%\\\/|\[\]()])/\\/g' | ## escape special characters @#$%\/|[]()
        awk -v chunk=${CHUNK} '{ print } NR % chunk == 0 { print "##" }' | ## split by "chunks"
        tr '\n' '|' | ## concat IDs with '|'
        sed 's/|##|/\n/g' | ## replace chunk separator with newline
        sed 's/^/^([^\\t]+\\t){8}(ID=|[^\\t]+;ID=)(/' | ## there must be 8 tab-separated columns before ID=XXX
        sed 's/|*$/)(;|$)/' > ${f_pattern} ## remove trailing '|', ID must end with either EOL or ';'
    ## check if no patterns. If no patterns, create empty output file
    local n_pattern=$(grep '(' ${f_pattern} | wc -l)
    if [ ${n_pattern} -eq 0 ]; then
        echo -n '' > ${fout}
    else
        printf '\n' >> ${f_pattern}
        while IFS= read -r line; do
            grep -nP "${line}" ${gff} >> ${fout}
        done < ${f_pattern}
    fi
    # rm ${f_pattern}
}

function subset_parent() {
    state "--subset_parent--"
    local gff fid fout gff_self dir_parent
    dir_parent=${dir_tmp}
    while (( "$#" )); do
        case "$1" in
            -g|--gff) gff="${2}"; shift;;
            -f|--fid) fid="${2}"; shift;;
            -o|--out|--fout) fout="${2}"; shift;;
            --gffself) gff_self="${2}"; shift;;
            -d|--dirparent) dir_parent="${2}"; shift;;
        esac
        shift
    done
    mkdir -p ${dir_parent}
    local f_pattern=$(mktemp -p ${dir_parent} parent.pattern.XXXXXX)
    ## generate gff for features in fid
    if [ -z ${gff_self} ]; then
        local tmp_gff_self=$(mktemp -p ${dir_parent} parent.self_gff.XXXXXX)
        subset_self --gff ${gff} --fid ${fid} --fout ${tmp_gff_self}
        gff_self=${tmp_gff_self}
    else
        ## unset, because tmp_gff_self may have been defined by earlier calls and we don't want to keep prev values
        unset tmp_gff_self
    fi
    cut -f9 ${gff_self} | ## extract attributes column
        grep -oP '(?<=^Parent=|;Parent=)[^;]+' | ## get PARENT attribute
        sort | uniq | ## remove duplicates
        perl -pe 's/(?=[\@#\$%\\\/|\[\]()])/\\/g' | ## escape special characters @#$%\/|[]()
        awk -v chunk=${CHUNK} '{ print } NR % chunk == 0 { print "##" }' | ## split by "chunks"
        tr '\n' '|' | ## concat IDs with '|'
        sed 's/|##|/\n/g' | ## replace chunk separator with newline
        sed 's/^/^([^\\t]+\\t){8}(ID=|[^\\t]+;ID=)(/' | ## there must be 8 tab-separated columns before ID=XXX
        sed 's/|*$/)(;|$)/' > ${f_pattern} ## remove trailing '|', ID must end with either EOL or ';'
    ## check if no patterns. If no patterns, create empty output file
    local n_pattern=$(grep '(' ${f_pattern} | wc -l)
    if [ ${n_pattern} -eq 0 ]; then
        echo -n '' > ${fout}
    else
        printf '\n' >> ${f_pattern}
        while IFS= read -r line; do
            grep -nP "${line}" ${gff} >> ${fout}
        done < ${f_pattern}
    fi
    # rm ${f_pattern} ${tmp_gff_self}
}


## in:
## /mnt/chaelab/rachelle/zzOtherzz/Zimin/pin3/phytozome14/20251110/462/v1/tmp2/subset_genome_around_blast_hits
## run:
## /mnt/chaelab/rachelle/src/subset_gff.sh --in /mnt/chaelab/rachelle/zzOtherzz/Zimin/tmp/gunzipped/462/Hvulgare_462_r1.gene_exons.gff3 --fid idA.txt --out tmp.gff --verbose --family --dir ./test_subset6

## TODO: use positional arguments for functions that are called by other functions
function subset_child() {
    state "---subset_child---"
    local gff fid fout dir_child
    dir_child=${dir_tmp}
    while (( "$#" )); do
        echo "1: ${1}; 2: ${2}"
        case "$1" in
            -g|--gff) gff="${2}"; shift;;
            -f|--fid) fid="${2}"; shift;;
            -o|--out|--fout) fout="${2}"; shift;;
            -d|--dirchild) dir_child="${2}"; shift;;
        esac
        shift
    done
    mkdir -p ${dir_child}
    echo "dir_child: ${dir_child}"
    local f_pattern=$(mktemp -p ${dir_child} child.pattern.XXXXXX)
    sort ${fid} | uniq | ## remove duplicates
        perl -pe 's/(?=[\@#\$%\\\/|\[\]()])/\\/g' | ## escape special characters @#$%\/|[]()
        awk -v chunk=${CHUNK} '{ print } NR % chunk == 0 { print "##" }' | ## split by "chunks"
        tr '\n' '|' | ## concat IDs with '|'
        sed 's/|##|/\n/g' | ## replace chunk separator with newline
        sed 's/^/^([^\\t]+\\t){8}(Parent=|[^\\t]+;Parent=)(/' | ## there must be 8 tab-separated columns before Parent=XXX
        sed 's/|*$/)(;|$)/' > ${f_pattern} ## remove trailing '|', ID must end with either EOL or ';'
    ## check if no patterns. If no patterns, create empty output file
    local n_pattern=$(grep '(' ${f_pattern} | wc -l)
    if [ ${n_pattern} -eq 0 ]; then
        echo -n '' > ${fout}
    else
        printf '\n' >> ${f_pattern}
        while IFS= read -r line; do
            grep -nP "${line}" ${gff} >> ${fout}
        done < ${f_pattern}
    fi
    # rm ${f_pattern}
}

# attribute_from_gff --gff ${gff_tmp} --fout ${fid_current} --attr ID
function attribute_from_gff() {
    state "------attribute_from_gff------"
    local fin="${1}"
    local fout="${2}"
    local attr="${3}"
    cut -f9 ${fin} | ## extract attributes column
        grep -oP "(?<=^${attr}=|;${attr}=)[^;]+" | ## get ID attribute
        sort | uniq > ${fout} ## remove duplicates
}

## generates txt file of <ID>\t<Parent ID> from GFF file
## ...let's just do this in python lmao
function parent_map() {
    state "------parent_map------"
    local gff fout
    while (( "$#" )); do
        case "$1" in
            -g|--gff) gff="${2}"; shift;;
            -o|-out|--fout) fout="${2}"; shift;;
        esac
        shift
    done
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); import regex as re; import data_manip as dm; attrs = [x.split('\t')[-1] for x in dm.splitlines('${gff}')]; match = [[re.search('(^ID=|;ID=)([^;]+)', x), re.search('(^Parent=|;Parent=)([^;]+)', x)] for x in attrs]; print(match); extracted = [['' if not x1 else x1.group(2), '' if not x2 else x2.group(2)] for x1,x2 in match]; dm.write_delim('${fout}', extracted, delim = '\t')"
    # awk 'BEGIN {FS="\t"; OFS="\t"} {id ~ /(?<=^ID=|;ID=)[^;]+/; parent ~ /(?<=^Parent=|;Parent=)[^;]+/; print id,parent}' ${gff} > ${fout}
}

function subset_ancestors() {
    state "---subset_ancestors---"
    local gff fid fout gff_self gff_seed
    while (( "$#" )); do
        case "$1" in
            -g|--gff) gff="${2}"; shift;;
            -f|--fid) fid="${2}"; shift;;
            -o|--out|--fout) fout="${2}"; shift;;
            --gffself) gff_self="${2}"; shift;;
        esac
        shift
    done
    ## generate gff for features in fid and set gff_seed
    local f_pattern=$(mktemp_f ancestors.pattern.XXXXXX)
    if [ -z ${gff_self} ]; then
        local tmp_gff_self=$(mktemp_f ancestors.self_gff.XXXXXX)
        subset_self --gff ${gff} --fid ${fid} --fout ${tmp_gff_self}
        gff_seed=${tmp_gff_self}
    else
        gff_seed=${gff_self}
    fi
    ## temporary paths
    local dir=${dir_tmp}/ancestors
    local fid_ancestors=${dir}/ancestors.fid.txt
    local fid_current=${dir}/current.fid.txt
    local fid_parent=${dir}/parent.fid.txt
    mkdir -p ${dir}
    ## get parents of parents etc. until no more parents are found
    while true; do
        gff_tmp=$(mktemp -p ${dir} gff_tmp.XXXXXX)
        # echo "gff_tmp: ${gff_tmp}"
        # echo "gff_seed: ${gff_seed}"
        ## get parents
        subset_parent --gff ${gff} --fid ${fid} --fout ${gff_tmp} \
                      --gffself ${gff_seed} -d ${dir}
        ## extract parent feature IDs
        attribute_from_gff ${gff_tmp} ${fid_current} ID
        ## append parent feature IDs to IDs of ancestors
        cat ${fid_current} >> ${fid_ancestors}
        ## check if current set of parent features has parents
        attribute_from_gff ${gff_tmp} ${fid_parent} Parent
        local n_parent=$(wc -l < ${fid_parent})
        ## exit if no further parents
        if [ ${n_parent} -eq 0 ]; then
            break
        fi
        ## update seed
        gff_seed=${gff_tmp}
    done
    ## get unique list of ancestor IDs
    local fid_ancestors_uniq=${dir}/ancestors.fid.uniq
    sort ${fid_ancestors} | uniq > ${fid_ancestors_uniq}
    ## subset input gff into fout
    subset_self --fid ${fid_ancestors_uniq} --fout ${fout} \
                --gff ${gff} -d ${dir}
    # ## remove temporary directory
    # rm -r ${dir}
}

function subset_descendants() {
    state "---subset_descendants---"
    local gff fid fout
    while (( "$#" )); do
        case "$1" in
            -g|--gff) gff="${2}"; shift;;
            -f|--fid) fid="${2}"; shift;;
            -o|--out|--fout) fout="${2}"; shift;;
        esac
        shift
    done
    local fid_seed=${fid}
    ## temporary paths
    local dir=${dir_tmp}/descendants
    local fid_descendant_parents=${dir}/descendant_parents.fid.txt
    local fid_current=${dir}/current.fid.txt
    local fid_new=${dir}/new.fid.txt
    local gff_tmp=${dir}/gff_tmp.gff
    mkdir -p ${dir}
    ## get children of children etc. until no more children are found
    cp ${fid_seed} ${fid_descendant_parents}
    echo "dir: ${dir}"
    while true; do
        ## get children
        subset_child --fid ${fid_seed} --fout ${gff_tmp} --gff ${gff} \
                     -d ${dir}
        ## check if there are any entries
        local n_children=$(wc -l < ${gff_tmp})
        ## exit if gff is empty
        if [ ${n_children} -eq 0 ]; then
            break
        fi
        ## extract child feature IDs
        attribute_from_gff ${gff_tmp} ${fid_current} ID
        ## sort current fids and check for new IDs
        sort -o ${fid_current} ${fid_current}
        comm -3 -1 ${fid_descendant_parents} ${fid_current} > ${fid_new}
        ## append new child feature IDs to IDs of descendants
        cat ${fid_new} >> ${fid_descendant_parents}
        ## update seed (this really only does anything at the first iteration)
        fid_seed=${fid_new}
    done
    ## get unique list of descendant IDs
    local fid_descendant_parents_uniq=${dir}/descendant_parents.fid.uniq.txt
    sort ${fid_descendant_parents} | uniq > ${fid_descendant_parents_uniq}
    ## subset input gff into fout
    ## ## this previously used subset_self, but it missed some descendants because they don't have IDs
    ## ## so this has been modified to use subset_child instead
    subset_child --gff ${gff} --fid ${fid_descendant_parents_uniq} \
                 --fout ${fout} -d ${dir}
    # ## remove temporary directory
    # rm -r ${dir}
}

files_to_concat=''

## retrieve self's entries
if [ ${PARENT} -ne 0 ] || [ ${SELF} -ne 0 ]; then
    state "Retrieving annotations for specified features"
    gff_self=${dir_tmp}/self.gff
    subset_self --fid ${FID} --fout ${gff_self} --gff ${GFF_IN}
    files_to_concat="${files_to_concat} ${gff_self}"
fi

## retrieve parent's entries
if [ ${PARENT} -ne 0 ]; then
    state "Retrieving annotations for parents"
    gff_parent=${dir_tmp}/parent.gff
    subset_parent --fid ${FID} --fout ${gff_parent} --gff ${GFF_IN} \
                  --gffself ${gff_self}
    files_to_concat="${files_to_concat} ${gff_parent}"
fi

## retrieve child's entries
if [ ${CHILD} -ne 0 ]; then
    state "Retrieving annotations for children"
    gff_child=${dir_tmp}/child.gff
    subset_child --fid ${FID} --fout ${gff_child} --gff ${GFF_IN}
    files_to_concat="${files_to_concat} ${gff_child}"
fi

## retrieve ancestor entries (recursively call subset_parent until no new entries are obtained)
if [ ${ANCESTORS} -ne 0 ]; then
    state "Retrieving annotations for ancestors"
    gff_ancestors=${dir_tmp}/ancestors.gff
    subset_ancestors --fid ${FID} --fout ${gff_ancestors} \
                     --gff ${GFF_IN} --gffself ${gff_self}
    files_to_concat="${files_to_concat} ${gff_ancestors}"
fi

## retrieve descendant entries (recursively call subset_child until no new entries are obtained)
if [ ${DESCENDANTS} -ne 0 ]; then
    state "Retrieving annotations for descendants"
    gff_descendants=${dir_tmp}/descendants.gff
    subset_descendants --fid ${FID} --fout ${gff_descendants} \
                       --gff ${GFF_IN}
    files_to_concat="${files_to_concat} ${gff_descendants}"
fi

echo "chkpt pre-concat"

## concat and sort output
if [ ${SORT} -eq 0 ]; then
    ## keep original entry order
    cat ${files_to_concat} | ## concatenate grepped output
        sed 's/:/\t/' | ## convert first ':' to a tab so we can parse grep line numbers as a column
        sort -k1,1n | ## sort by original line numbers
        uniq | ## de-duplicate any repeated entries
        cut -f2- > ${GFF_OUT} ## remove grep line numbers
else
    ## sort by chromosome, start, and end coordinates
    cat ${files_to_concat} | ## concatenate grepped output
        sed 's/:/\t/' | ## convert first ':' to a tab so we can parse grep line numbers as a column
        cut -f2- | ## remove grep line numbers
        sort -k1,1 -k4,4n -k5,5n | ## sort by chromosome, start, end
        uniq > ${GFF_OUT} ## de-duplicate any repeated entries
fi

echo "chkpt end"

exit 0

# ## compare files
# if diff file1 file2 > /dev/null
# then
#     echo "No difference"
# else
#     echo "Difference"
# fi
