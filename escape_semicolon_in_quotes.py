#!/usr/bin/python3

import sys
import fileinput
## fileinput allows this script to take either a file or stdin

## this is only intended for use with non-overlapping quotes!!!
## prints to stdout

try:
    for line in fileinput.input(sys.argv[1:]):
        str_out = ''
        open_quote = False
        for c in line:
            if open_quote:
                if c == '"': open_quote = False
                if c == ';':
                    str_out += '%3B'
                    continue
            else:
                if c == '"': open_quote = True
            str_out += c
        print(str_out.rstrip('\n\r'))
except IOError:
    ## stdout is closed, no point in continuing
    ## attempt to close explicitly to prevent cleanup problems
    sys.stdout.close()
    sys.stderr.close()
