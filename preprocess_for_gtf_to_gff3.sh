#!/bin/bash

fin=$1

sed -r 's/([^;])$/\1;/' ${fin} | /mnt/chaelab/rachelle/src/escape_semicolon_in_quotes.py
