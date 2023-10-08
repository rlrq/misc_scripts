#!/bin/bash

## prints comma-separate output (likely going to be a filename, but doesn't have to be)
make_all(){
    ## arg1: function that takes 1 input and prints an output
    ## args2 onward: elements to pass to arg1 function
    local f=${1}
    local tmp=( ${@} )
    local arr=${tmp[@]:1}
    local outpt=""
    for e in ${arr[@]}; do
        outpt="${outpt},$(${f} ${e})"
    done
    outpt=$(echo ${outpt} | sed "s/^,//")
    echo ${outpt}
}

## make_all function, but all elements are additionally enclosed in single quotes
make_all_quoted(){
    ## same arguments as make_all
    local outpt=$(make_all ${@})
    outpt=$(echo ${outpt} | sed "s/,/','/g")
    echo "'${outpt}'"
}
