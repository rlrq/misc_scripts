#!/bin/bash

## although mafft is referred to as ALN, this script will ONLY work for mafft (cuz we're using the --add parameter to chaining alignments) lol. ALN is mostly there in case there are different versions of these tools
DIR_DEFAULT=$(pwd)
ALN_DEFAULT=/usr/bin/mafft
TRANSLATE_DEFAULT=/mnt/chaelab/rachelle/src/translate.py
PEP2CDS_DEFAULT=/mnt/chaelab/rachelle/src/aln_pep2cds.py
PREFIX_DEFAULT=pep2cds2geneAln
KEEP_CDS_DEFAULT=0

## if FOUT_CDS and/or FOUT_PEP are not provided, temporary files will be used to store the alignments and will be removed upon completion of programme
while (( "$#" )); do
    case "$1" in
        --dir) DIR="${2}";; ## output directory
        --pref|--prefix) PREFIX="${2}";; ## output prefix, used if --aln-out and/or --tree-out are not specified
        --cds) FA_CDS="${2}";; ## fasta of CDS (coding strand)
        --gene) FA_GENE="${2}";; ## fasta of gene (either strand is fine)
        --aln-pep) ALN_ARGS_PEP="${2}";; ## arguments for alignment (peptide)
        # --aln-cds) ALN_ARGS_CDS="${2}";; ## arguments for alignment (cds)
        --aln-gene) ALN_ARGS_GENE="${2}";; ## arguments for alignment (gene/genome)
        --mafft) ALN_ARGS="${2}";; ## arguments for alignment
        -o|--out|--out-gene) FOUT_GENE="${2}";; ## output genomic alignment file
        --out-cds) FOUT_CDS="${2}";; ## output CDS alignment file
        --out-pep) FOUT_PEP="${2}";; ## output peptide alignment file
        --which-aln) ALN="${2}";; ## path to alignment binary
        --which-translate) TRANSLATE="${2}";; ## path to translation script (takes 2 arguments: input and output paths)
        --which-pep2cds) PEP2CDS="${2}";; ## path to pep2cds script (takes 3 arguments: peptide alignment, fasta of CDS [should have identical sequence IDs to pep alignment], and output path)
        --keep-cds) KEEP_CDS=1;; ## keep CDS sequences in gene/genomic alignment (if not raised, sequence IDs MUST NOT BE THE SAME between FA_CDS and FA_GENE)
    esac
    shift
done

## set some variables
DIR="${DIR:-${DIR_DEFAULT}}"
PREFIX="${PREFIX:-${PREFIX_DEFAULT}}"
ALN="${ALN:-${ALN_DEFAULT}}"
TRANSLATE="${TRANSLATE:-${TRANSLATE_DEFAULT}}"
PEP2CDS="${PEP2CDS:-${PEP2CDS_DEFAULT}}"
KEEP_CDS="${KEEP_CDS:-${KEEP_CDS_DEFAULT}}"
if [ $(basename ${ALN}) == "mafft" ]; then
    aln_suff='mafft'
else
    aln_suff='aln'
fi
FOUT_GENE="${FOUT_GENE:-${DIR}/${PREFIX}.gene.${aln_suff}.fasta}"
# FOUT_CDS="${FOUT_CDS:-${DIR}/${PREFIX}.CDS.${aln_suff}.fasta}"
# FOUT_PEP="${FOUT_OUT:-${DIR}/${PREFIX}.pep.${aln_suff}.fasta}"

## print some stuff for checking
echo "Gene/genome alignment output: ${FOUT_GENE}"
echo "CDS alignment output: ${FOUT_CDS}"
echo "Peptide alignment output: ${FOUT_PEP}"

## create output directory
mkdir -p ${DIR}

## create temporary directory
TMPDIR=$(mktemp -d --tmpdir=${DIR})
mkdir -p ${TMPDIR}
trap 'rm -rf -- "${TMPDIR}"' 0 1 2 15 EXIT ## delete temporary directory on exit or interruption

## translate
tmp_pep=$(mktemp -p ${TMPDIR} "${PREFIX}.XXX")
${TRANSLATE} ${FA_CDS} ${tmp_pep}

## align peptides
tmp_pep_aln=$(mktemp -p ${TMPDIR} "${PREFIX}.XXX")
${ALN} ${ALN_ARGS_PEP} ${tmp_pep} > ${tmp_pep_aln}
## write to user-provided path
if [[ -n ${FOUT_PEP} ]]; then
    mv ${tmp_pep_aln} ${FOUT_PEP}
else
    ## remove unaligned pep sequences
    rm ${tmp_pep}
fi

## convert peptide alignment to CDS alignment
tmp_cds_aln=$(mktemp -p ${TMPDIR} "${PREFIX}.XXX")
${PEP2CDS} ${tmp_pep_aln} ${FA_CDS} ${tmp_cds_aln}
## write to user-provided path
if [[ -n ${FOUT_CDS} ]]; then
    mv ${tmp_cds_aln} ${FOUT_CDS}
else
    ## remove pep alignment
    rm ${tmp_pep_aln}
fi

## align genomic sequences with CDS
tmp_gene_aln=$(mktemp -p ${TMPDIR} "${PREFIX}.XXX")
${ALN} ${ALN_ARGS_GENE} --add ${FA_GENE} ${tmp_cds_aln} > ${tmp_gene_aln}
## keep/remove CDS sequences
if [ ${KEEP_CDS} -eq 1 ]; then
    mv ${tmp_gene_aln} ${FOUT_GENE}
else
    python3 -c "from Bio import SeqIO; cds_seqids = set(record.description for record in SeqIO.parse('${FA_CDS}', 'fasta')); aln_to_keep = [record for record in SeqIO.parse('${tmp_gene_aln}', 'fasta') if record.description not in cds_seqids]; SeqIO.write(aln_to_keep, '${FOUT_GENE}', 'fasta')"
    rm ${tmp_gene_aln}
fi
