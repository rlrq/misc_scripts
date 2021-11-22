#!/bin/bash

## tested with sha256sum

if [ -z ${3} ]; then
    echo "--missing argument(s)--"
    exit 1
fi

sum=$(${1} ${2} | tr ' ' '\t' | cut -f1)
check=$(printf "${sum}\n${3}" | uniq | wc -l)

if [[ ${check} -eq 1 ]]; then
    echo "--verified--"
else
    echo "--mismatch--"
fi
