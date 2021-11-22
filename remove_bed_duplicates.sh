#!/bin/bash

# removes duplicates in bed file and resorts it by chr coordinates
for arg; do
    echo ${arg}
    sort -u "${arg}" | bedtools sort -i - > ${arg}_temp.bed
    mv ${arg}_temp.bed ${arg}
done
