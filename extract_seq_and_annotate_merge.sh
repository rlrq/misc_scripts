#!/bin/bash

# Subsets GFF3 annotation by feature IDs, includes subfeatures.
# Extracts relevant sequences from larger FASTA file.

# Takes a file of feature IDs.
# Extracts feature ranges from GFF3 file.
# Extends feature range by 'extend' bp and merges all overlapping or bookended sequences to obtain extended range.
# Intersects extended range with GFF3 file again to obtain extended features.
# If any feature extends beyond the extended range, futher extend the range to include the entire feature.
# Extends further extended range by 'flank' bp to obtain flanked range.
# Extracts flanked range sequences from FASTA file, named by {molecule}:{original start (1-index)}.{original end (inclusive)}.
# Subtracts the lower value of the relevant flanked range from all coordinates in extended features.
# Updates the molecule of all extended features to that of the relevant flanked range.
# Writes GFF3 subsetted from original and updated with modified coordinates to file.
# writes FASTA of flanked range sequences to file.
# Writes GFF3 subsetted from original to file (coordinates per original GFF3 file) if gff_original is provided.

EXTEND_BP_DEFAULT=10000
FLANK_BP_DEFAULT=500
EXTEND_PARENT_DEFAULT=1
MEMSAVE_DEFAULT=1
SORT_DEFAULT=0
THRESHOLD_DEFAULT=511999999 ## maximum chromosome size for bedtools without -sorted option is 512 Mb (https://bedtools.readthedocs.io/en/latest/)

while (( "$#" )); do
    case "$1" in
        --gff) GFF="$(realpath ${2})"; shift;; ## input GFF3 file
        --fasta) FASTA="$(realpath ${2})"; shift;; ## input FASTA file
        --id) F_FID="$(realpath ${2})"; USE_ID=1; USE_RANGE=0; shift;; ## input TXT file of feature IDs (1 ID per line)
        --range) F_RANGE="$(realpath ${2})"; USE_ID=0; USE_RANGE=1; shift;; ## input BED file of ranges (1 per line; format: <chrom>\t<start>\t<end>)
        --extend) EXTEND_BP="${2}"; shift;; ## number of bp to extend for feature subsetting
        --flank) FLANK_BP="${2}"; shift;; ## number of bp to extend beyond extended features for FASTA subsetting
        --gff-out) GFF_OUT="$(realpath ${2})"; shift;; ## output GFF3 file
        --fasta-out) FASTA_OUT="$(realpath ${2})"; shift;; ## output FASTA file
        ## output BED file containing mapping coords to sequence IDs in FASTA_OUT
        --bed-out) BED_OUT="$(realpath ${2})"; shift;;
        ## write GFF3 subset with original coords to this file
        --gff-original|--gff-og) GFF_ORIGINAL="$(realpath ${2})"; shift;;
        ## move temporary directory to DIR_TMP_OUT upon exit instead of deleting it
        --tmp-out) DIR_TMP_OUT="$(realpath ${2})"; shift;;
        ## extend around coordinates of highest parent for each feature, not the feature itself
        --parent) EXTEND_PARENT=1;;
        --self) EXTEND_PARENT=0;; ## extend around coordinates of feature, not its highest parent
        --memsave) MEMSAVE=1;; ## use memsave mode
        --no-memsave) MEMSAVE=0;; ## don't use memsave mode
        ## sorts GFF3 file and then pass -sorted argument and sorted GFF3 file to bedtools intersect
        --sort) SORT=1;;
        ## chromosome size threshold for automatically setting SORT=1
        --sort-threshold|--threshold) THRESHOLD="${2}"; shift;;
    esac
    shift
done

## set some variables with defaults
FLANK_BP="${FLANK_BP:-${FLANK_BP_DEFAULT}}"
EXTEND_BP="${EXTEND_BP:-${EXTEND_BP_DEFAULT}}"
EXTEND_PARENT=${EXTEND_PARENT:-${EXTEND_PARENT_DEFAULT}}
MEMSAVE=${MEMSAVE:-${MEMSAVE_DEFAULT}}
SORT=${SORT:-${SORT_DEFAULT}}
THRESHOLD=${THRESHOLD:-${THRESHOLD_DEFAULT}}

## check that required arguments are provided
if [ -z ${GFF} ]; then
    echo "Input GFF3 file is required."
    exit 1
elif [ -z ${FASTA} ]; then
    echo "Input FASTA file is required."
    exit 1
elif [ -z ${F_FID} ] && [ -z ${F_RANGE} ]; then
    echo "Input feature ID file (1 ID per line) or coordinate range file (1 range per line) is required."
    exit 1
elif [ -z ${GFF_OUT} ]; then
    echo "Output GFF file is required."
    exit 1
elif [ -z ${FASTA_OUT} ]; then
    echo "Output FASTA file is required."
    exit 1
fi

## set temporary directory
dir_tmp=$(mktemp -d)
function mktemp_f { echo "$(mktemp -p ${dir_tmp})"; }

## cleanup
function cleanup_delete { echo "Removing temporary directory '${dir_tmp}'"; rm -r ${dir_tmp}; }
function cleanup_move { echo "Moving temporary directory '${dir_tmp}' to '${DIR_TMP_OUT}'"; mv ${dir_tmp} ${DIR_TMP_OUT}; }
if [ -z ${DIR_TMP_OUT} ]; then
    trap cleanup_delete EXIT
    trap cleanup_delete SIGINT
else
    trap cleanup_move EXIT
    trap cleanup_move SIGINT
fi

## set other variables
if [ ${MEMSAVE} -eq 0 ]; then
    memsave_arg=''
else
    memsave_arg='--memsave'
fi
# if [ ${EXTEND_PARENT} -eq 0 ]; then
#     parent_arg=''
# else
#     parent_arg='--parent'
# fi

## set SORT=1 if largest chromosome size exceeds THRESHOLD
largest_chom_size=$(cut -f5 ${GFF} | sort -k1,1n | tail -1)
if [ ${largest_chom_size} -ge ${THRESHOLD} ]; then
    SORT=1
fi


##########################################################
echo "---------------------------------------------------"
echo "####        dir_tmp: ${dir_tmp}       ####"
echo "---------------------------------------------------"
##          Retrieve annotations for features           ##
##########################################################

if [ ${USE_ID} -eq 1 ]; then
    echo "F_FID: ${F_FID}"
    ## get gff annotations for relevant features
    tmp_gff=$(mktemp_f)
    if [ ${EXTEND_PARENT} -eq 1 ]; then
        ## if using parents' ranges
        echo "---Retrieving coordinates of features and parents---"
        /mnt/chaelab/rachelle/src/subset_gff.sh --in ${GFF} --fid ${F_FID} --out ${tmp_gff} --quiet --self --parent
        # /mnt/chaelab/rachelle/src/extract_gff_features.py ${GFF} ${F_FID} ${tmp_gff} --fmt GFF --parent ${memsave_arg}
    else
        ## if using ranges of provided features only
        echo "---Retrieving coordinates of features---"
        /mnt/chaelab/rachelle/src/subset_gff.sh --in ${GFF} --fid ${F_FID} --out ${tmp_gff} --quiet --self
        # /mnt/chaelab/rachelle/src/extract_gff_features.py ${GFF} ${F_FID} ${tmp_gff} --fmt GFF --self ${memsave_arg}
    fi
fi

##########################################################
echo "---Subsetting GFF for features in extended range---"
##########################################################

## get extended ranges by extending using EXTEND_BP
tmp_bedext=$(mktemp_f)
if [ ${USE_ID} -eq 1 ]; then
    # awk -v e=${EXTEND_BP} 'BEGIN {FS="\t"; OFS="\t"} {if($4-1-e<0)$4=e; print $1,$4-1-e,$5+e}' ${tmp_gff} |
    #     sort -k1,1 -k2,2n | bedtools merge -i stdin > ${tmp_bedext}
    awk -v e=${EXTEND_BP} 'BEGIN {FS="\t"; OFS="\t"} {if($4-1-e<0)$4=e+1; print $1,$4-1-e,$5+e}' ${tmp_gff} |
        sort -k1,1 -k2,2n | bedtools merge -i stdin > ${tmp_bedext}
elif [ ${USE_RANGE} -eq 1 ]; then
    ## create extended range from input bed file if --range
    awk -v e=${EXTEND_BP} 'BEGIN {FS="\t"; OFS="\t"} {if($2-e<0)$2=e; print $1,$2-e,$3+e}' ${F_RANGE} |
        sort -k1,1 -k2,2n | bedtools merge -i stdin > ${tmp_bedext}
fi

## created sorted GFF if --sort
if [ ${SORT} -eq 1 ]; then
    tmp_gff_sorted=$(mktemp_f)
    grep -v '^#' ${GFF} | sort -k1,1 -k4,4n -k5,5n > ${tmp_gff_sorted}
    GFF=${tmp_gff_sorted} ## replace $GFF with a temporary file
fi

## intersect extended range with GFF
tmp_gffext=$(mktemp_f)
if [ ${SORT} -eq 0 ]; then
    bedtools intersect -a ${tmp_bedext} -b ${GFF} -wb | awk '$3!="chromosome"' > ${tmp_gffext}
else
    bedtools intersect -a ${tmp_bedext} -b ${GFF} -wb -sorted | awk '$3!="chromosome"' > ${tmp_gffext}
fi

## subset GFF for final set of features
tmp_fid=$(mktemp_f)
gff_og=$(mktemp_f)
## only keep IDs that are top-level (i.e. do not have parents)
cut -f12 ${tmp_gffext} | grep -vP '^Parent=|;Parent=' | grep -oP '^ID=[^;]+|;ID=[^;]+' | sed 's/^ID=//' | sort | uniq > ${tmp_fid}
/mnt/chaelab/rachelle/src/subset_gff.sh --in ${GFF} --fid ${tmp_fid} --out ${gff_og} --quiet --family
# /mnt/chaelab/rachelle/src/extract_gff_features.py ${GFF} ${tmp_fid} ${gff_og} --fmt GFF --include-cousins ${memsave_arg}
if [ -n ${GFF_ORIGINAL} ]; then
    mv ${gff_og} ${GFF_ORIGINAL}
    gff_og=${GFF_ORIGINAL}
fi


##############################################################
echo "---Extracting sequence for extended range with flank---"
##############################################################

## append extracted feature ranges to tmp_bedext
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$7-1,$8}' ${tmp_gffext} >> ${tmp_bedext}
## further extend extended range by merging to fully include any feature that partially overlaps with the original extended range
## and then extend again using FLANK_BP
## also, assign new "molecule" ID with pattern "<original_molecule>:<start (1-index inclusive)>..<end (inclusive)>"
tmp_bedext2=$(mktemp_f)
sort -k1,1 -k2,2n ${tmp_bedext} | bedtools merge -i stdin |
    awk -v f=${FLANK_BP} 'BEGIN {FS="\t"; OFS="\t"} {if($2-f<0)$2=f; print $1,$2-f,$3+f,$1":"$2+1-f".."$3+f}' > ${tmp_bedext2}
if [ -n ${BED_OUT} ]; then
    mv ${tmp_bedext2} ${BED_OUT}
    tmp_bedext2=${BED_OUT}
fi

## extract FASTA sequence
python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); import index_fasta as ifasta; import data_manip as dm; import fasta_manip as fm; fa = ifasta.IndexedFasta('${FASTA}'); ranges = [x.split('\t') for x in dm.splitlines('${tmp_bedext2}')]; ranges = [(x[0],int(x[1]),int(x[2]),x[3]) for x in ranges]; seqs = {new_id: fa[mol][start:end] for mol,start,end,new_id in ranges}; fm.dict_to_fasta(seqs, '${FASTA_OUT}')"


###################################################################
echo "---Generating GFF with adjusted ranges in subsetted FASTA---"
###################################################################

## update GFF annotation with adjusted ranges in subsetted FASTA_OUT
## first intersect gff_og with tmp_bedext2 to match gff entries to the merged ranges
tmp_bed_gff=$(mktemp_f)
if [ ${SORT} -eq 0 ]; then
    bedtools intersect -a ${tmp_bedext2} -b ${gff_og} -loj |
        awk -v f=${FLANK_BP} 'BEGIN {FS="\t"; OFS="\t"} {print $4,$6,$7,$8-$2,$9-$2,$10,$11,$12,$13}' > ${GFF_OUT}

else
    bedtools intersect -a ${tmp_bedext2} -b ${gff_og} -loj -sorted |
        awk -v f=${FLANK_BP} 'BEGIN {FS="\t"; OFS="\t"} {print $4,$6,$7,$8-$2,$9-$2,$10,$11,$12,$13}' > ${GFF_OUT}
fi




##########################################################
echo "---------------------------------------------------"
##                         EXIT                         ##
##########################################################

# return

