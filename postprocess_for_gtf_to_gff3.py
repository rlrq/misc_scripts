#!/usr/bin/python3

import re
import sys
import fileinput
## fileinput allows this script to either take a file or stdin

## this script replaces default genometools gtf_to_gff3 ID/Parent tag values with relevant gene_id or transcript_id values
## afaik only 2 types of features get ID tags: gene and mRNA, but we can extend this script if we find more
## and 2 types of features (exon and mRNA) get Parent tags

def get_gene_id_value(s):
    return re.search("(?<=;gene_id=)[^;]+|(?<=^gene_id=)[^;]+", s)[0]
def get_transcript_id_value(s):
    return re.search("(?<=;transcript_id=)[^;]+|(?<=^transcript_id=)[^;]+", s)[0]

ID_pattern = re.compile(r"[^;]ID=([^;]+)")
Parent_pattern = re.compile(r"[^;]Parent=([^;]+)")

def sub_ID(s, val):
    pre_ID = 'ID=' if re.match("ID=", s) else re.search("^.*?;ID=", s)[0]
    post_ID = re.search('^' + pre_ID + "([^;]+)(.*?$)", s)[2]
    return pre_ID + val + post_ID
def sub_Parent(s, val):
    pre_Parent = 'Parent=' if re.match("Parent=", s) else re.search("^.*?;Parent=", s)[0]
    post_Parent = re.search('^' + pre_Parent + "([^;]+)(.*?$)", s)[2]
    return pre_Parent + val + post_Parent

try:
    for line in fileinput.input(sys.argv[1:]):
        ## just print comment lines as-is
        if line[0] == '#':
            print(line.rstrip('\n\r'))
            continue
        ## process non-comment lines
        split_line = line.rstrip('\n\r').split('\t')
        feature = split_line[2]
        attributes = split_line[-1]
        if feature == "gene":
            attributes = sub_ID(attributes, get_gene_id_value(attributes))
        elif feature == "mRNA":
            attributes = sub_ID(attributes, get_transcript_id_value(attributes))
            attributes = sub_Parent(attributes, get_gene_id_value(attributes))
        elif feature == "exon":
            attributes = sub_Parent(attributes, get_transcript_id_value(attributes))
        split_line[-1] = attributes
        print('\t'.join(split_line))
except IOError:
    ## stdout is closed, no point in continuing
    ## attempt to close explicitly to prevent cleanup problems
    sys.stdout.close()
    sys.stderr.close()
