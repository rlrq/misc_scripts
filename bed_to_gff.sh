#!/bin/bash

bed_in=$1
gff_out=$2

# intersect
# nlr_f=/mnt/chaelab/rachelle/NLR/athaliana_nblrr_20kbmerge.bed
# intersect_f=/mnt/chaelab/rachelle/tDNA/bed/nblrr_bg_gff.bed
# bedtools intersect -a $bed_f -b $nlr_f -wa > $intersect_f
# intersect_gff=/mnt/chaelab/rachelle/tDNA/bed/nblrr_bg_gff.gff

# convert back to original gff format
awk '{print $1 "\t" $4 "\t" $5 "\t" $2 " \t" $3 "\t" $6 " \t" $7 "\t" $8 "\t" $9}' $bed_in > $gff_out
