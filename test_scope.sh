#!/bin/bash

while (( "$#" )); do
    case "$1" in
        -a) A="${2}"; shift;;
        -b) B="${2}"; shift;;
    esac
    shift
done

echo "${A}"
echo "${B}"

function f1() {
    echo "f1 1: ${1}"
}

function f2() {
    echo "f2 1: ${1}"
    f1 "f2_to_f1"
}

function f3() {
    local A keyword_a
    while (( "$#" )); do
        case "$1" in
            -i) A="${2}"; shift;;
            -a) keyword_a="${2}"; shift;;
        esac
        shift
    done
    echo "f3 A (-i): ${A}"
    echo "f3 keyword_a (-a): ${keyword_a}"
}

function f4() {
    echo "f4 1: ${1}"
    f3 -i "f4_to_f3_i" -a "f4_to_f3_a"
}

function f5() {
    local A B
    while (( "$#" )); do
        case "$1" in
            -a) A="${2}"; shift;;
            -b) B="${2}"; shift;;
        esac
        shift
    done
    echo "f5 A: ${A}"
    echo "f5 B: ${B}"
}

function f6() {
    local A B
    while (( "$#" )); do
        case "$1" in
            -a) A="${2}"; shift;;
            -b) B="${2}"; shift;;
        esac
        shift
    done
    echo "f6 A: ${A}"
    echo "f6 B: ${B}"
    f5 -a "f6_to_f5_a" -b "f6_to_f5_b"
}

function f7() {
    f5 -a "f7_to_f5_a" -b "f7_to_f5_b"
}

echo "--testing f1--"
f1 "f1_arg1"
echo "--testing f2--"
f2 "f2_arg1"
echo "--testing f3--"
f3 -i "f3_arg_i" -a "f3_arg_a"
echo "--testing f4--"
f4 "f4_arg1"
echo "--testing f5--"
f5 -a "f5_a" -b "f5_b"
echo "global A: ${A}"
echo "global B: ${B}"
echo "--testing f6--"
f6 -a "f6_a" -b "f6_b"
echo "--testing f7--"
f7
