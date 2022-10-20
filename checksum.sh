#!/bin/bash

## tested with sha256sum

if [ -z ${3} ]; then
    echo "--missing argument(s)--"
    echo "usage: checksum <algorithm> <fname> <sum>"
    echo "(e.g. 'checksum md5sum interproscan-5.56-89.0-64-bit.tar.gz f3288327a6c7be1ff42220c3b7a1de43')"
    exit 1
fi

sum=$(${1} ${2} | tr ' ' '\t' | cut -f1)
check=$(printf "${sum}\n${3}" | uniq | wc -l)

if [[ ${check} -eq 1 ]]; then
    echo "--verified--"
else
    echo "--mismatch--"
fi
