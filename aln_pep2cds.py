#!/usr/bin/python3

import re
import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from fasta_manip import map_cds_to_pep_aln

aln_pep, fa_cds, aln_fout = sys.argv[1:4]
pep_pattern, cds_pattern = (".+", ".+") if len(sys.argv) < 5 else sys.argv[4:6]

map_cds_to_pep_aln(aln_fout, aln_pep, fa_cds,
                   name_map = lambda pep_id, cds_id: (re.search(pep_pattern, pep_id).group(0) == \
                                                      re.search(cds_pattern, cds_id).group(0)))
