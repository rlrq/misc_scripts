get_column_by_name () {
    ## outputs to stdout
    local fname=${1}
    local colname=${2}
    awk -v col=${colname} -v FS='\t' 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR>1{print $c}' ${fname}
}
