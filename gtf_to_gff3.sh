#!/bin/bash

## wraps a couple of commands to include preprocessing before conversion to gff3

fin=$1

/mnt/chaelab/rachelle/src/preprocess_for_gtf_to_gff3.sh ${fin} |
    /mnt/chaelab/rachelle/programmes/GenomeTools/genometools-1.6.2/bin/gt gtf_to_gff3 ${@:2} |
    /mnt/chaelab/rachelle/src/postprocess_for_gtf_to_gff3.py
