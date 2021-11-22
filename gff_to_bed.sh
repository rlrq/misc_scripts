#!/bin/bash

gff_in=$1
bed_out=$2
func=$3

# gff_in=/mnt/chaelab/shared/genomes/a_thaliana/fasta/TAIR10/features/TAIR10_GFF3_genes.gff
# bed_out=/mnt/chaelab/rachelle/tDNA/bed/gff_bed.bed

if [ -z "${func}" ]; then
    echo "leaving chromosomes as is"
    paste -d"\t" <(awk -F '\t' '{print $1 "\t" $4-1 "\t" $5}' ${gff_in}) <(cut -d$'\t' -f 9 ${gff_in} | sed 's/.*ID=//g' | sed 's/;.*$//g' | sed 's/^.*=.*$/./g') <(awk -F '\t' '{print $6 "\t" $7 "\t" $2 " \t" $3 "\t" $8 "\t" $9}' $gff_in) > $bed_out
elif [ $func == 1 ]; then
    echo "removing chromosome"
    paste -d"\t" <(awk -F '\t' '{print $1 "\t" $4-1 "\t" $5}' ${gff_in}) <(cut -d$'\t' -f 9 ${gff_in} | sed 's/.*ID=//g' | sed 's/;.*$//g' | sed 's/^.*=.*$/./g') <(awk -F '\t' '{print $6 "\t" $7 "\t" $2 " \t" $3 "\t" $8 "\t" $9}' $gff_in) | awk -F "\t" '{gsub("Chr","chr",$1)}1' OFS="\t" | awk -F "\t" '$5 != "chromosome" {print}' OFS="\t" > $bed_out
elif [ $func == 2 ]; then
    echo "retaining chromosome"
    # retain chromosome
    paste -d"\t" <(awk -F '\t' '{print $1 "\t" $4-1 "\t" $5}' ${gff_in}) <(cut -d$'\t' -f 9 ${gff_in} | sed 's/.*ID=//g' | sed 's/;.*$//g' | sed 's/^.*=.*$/./g') <(awk -F '\t' '{print $6 "\t" $7 "\t" $2 " \t" $3 "\t" $8 "\t" $9}' $gff_in) | awk -F "\t" '{gsub("Chr","chr",$1)}1' OFS="\t" > $bed_out
fi
