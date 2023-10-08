#!/bin/bash

function getseq_to_bed {
    sed 's/|/\t/g' $1 | cut -f1,5 | sed 's/-/\t/g' | sed 's/^_R_//g' | paste -d '\t' - $1 > $2
}

# function get_vdw_gene_overlapping_bed {
#     local vdw_contigs="/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"
#     local vdw_gff="/mnt/chaelab/rachelle/data/anna_lena/vdw_gff3.bed"
#     awk '$8=="gene"' ${vdw_gff} | bedtools intersect -wao -a $1 -b stdin | awk '{if ($8!=".") print $8}' > $2
#     /mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh -g $2 -f gene -b ${vdw_gff} -a ref --fasta ${vdw_contigs} --out $3 -r ref --adjust-dir --no-bed
# }

# ## args:
# ##   <bed of query for overlap> <output file> <nonref gff_bed> <nonref fasta>
# function get_nonref_gene_overlapping_bed {
#     /mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh -g $( awk '$8=="gene"' $3 | bedtools intersect -wao -a $1 -b stdin | awk '{if ($8!=".") print $8}' | tr '\n' ',' ) -f gene -b ${vdw_gff} -a ref --fasta $4 --out $2 -r ref --adjust-dir --no-bed
# }

## args:
##   <bed of query for overlap> <output file> <nonref gff_bed> <nonref fasta> <feature>
function get_nonref_feature_overlapping_bed {
    feature="${5:-gene}"
    /mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh -g $( awk '$8=="gene"' $3 | bedtools intersect -wao -a $1 -b stdin | awk '{if ($8!=".") print $8}' | tr '\n' ',' ) -f ${feature} -b $3 -a ref --out $2 -r $4 --adjust-dir --no-bed
}

## args:
##   1. <file of seqnames (output from getSeq)>
##   2. <file to write bed of seqnames>
##   3. <output file>
##   4. <nonref gff_bed>
##   5. <nonref fasta>
##   6. <feature>
function get_nonref_feature_from_name {
    getseq_to_bed $1 $2
    get_nonref_feature_overlapping_bed $2 $3 $4 $5 $6
}

function get_nonref_gene_from_name {
    getseq_to_bed $1 $2
    get_nonref_feature_overlapping_bed $2 $3 $4 $5 gene
}

function get_vdw_feature_overlapping_bed {
    get_nonref_feature_overlapping_bed $1 $2 "/mnt/chaelab/rachelle/data/anna_lena/vdw_gff3.bed" "/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta" $3
}

function get_vdw_gene_overlapping_bed {
    get_nonref_feature_overlapping_bed $1 $2 "/mnt/chaelab/rachelle/data/anna_lena/vdw_gff3.bed" "/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta" gene
}

function get_vdw_feature_from_name {
    getseq_to_bed $1 $2
    get_vdw_feature_overlapping_bed $2 $3 $4
}

function get_vdw_gene_from_name {
    getseq_to_bed $1 $2
    get_vdw_feature_overlapping_bed $2 $3 gene
}

function get_seq_from_gid {
    local fname=${1}
    local fout=${2}
    local gids=${3}
    local feature_type=${4}
    python3 -c "import re; from Bio import SeqIO; ids = set('${gids}'.split(',')); records = [record for record in SeqIO.parse('${fname}', 'fasta') if re.search('(?<=\|).+?(?=\|${feature_type}\|)', record.id).group(0) in ids]; SeqIO.write(records, '${fout}', 'fasta')"
}
