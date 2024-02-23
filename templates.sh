#!/bin/bash

## parallel process
## from: https://unix.stackexchange.com/a/436713
N=4
sample_array=( 1 2 3 )
for i in ${sample_array[@]}; do
    (
        ## stuff to parallelise in here
        echo "dummy"
    ) &
    ## allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge ${N} ]]; then
        ## now there are $N jobs already runnin, so wait here for any job
        ## to be finished so there is a place to start next one
        wait -n
    fi
done
## no more jobs to be started but wait for pending jobs
## (all need to be finished)
wait
echo "finished"

## read file
inpt_fname=/path/to/file
while read line; do
    ## process stuff here (i think ${line} doesn't include newline)
    echo "${line}"
done <${inpt_fname}

## read file w/o interpreting backslash (-r) and w/o trimming leading/trailing whitespace
inpt_fname=/path/to/file
while IFS= read -r line; do
    ## process stuff here (i think ${line} doesn't include newline)
    echo "${line}"
done <${inpt_fname}

## more fun stuff on reading into loops: https://www.cyberciti.biz/faq/unix-howto-read-line-by-line-from-file/
