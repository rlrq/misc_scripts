#!/bin/bash

params=("$@")
n="$1"
command=("${params[@]:1}")

for i in $(seq ${n}); do
    eval "${command[@]}"
done

