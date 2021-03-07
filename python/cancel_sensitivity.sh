#! /bin/bash
[[ $# -ne 1 ]] && { echo Usage: $0 BATCH_NAME; exit 1; }
BATCH_NAME=$1
nb=${#BATCH_NAME}

nums=($(squeue -o "%.8i" -u $(whoami)))
names=($(squeue -o "%.40j" -u $(whoami)))
[[ ${#nums[@]} -ne ${#names[@]} ]] && exit 1

njobs=${#nums[@]}
njobs=$((njobs-1))
for i in `seq 2 $njobs`
do
    i_=$((i-1))
    num=${nums[$i_]}
    name=${names[$i_]}
    if [ ${name:0:$nb} == $BATCH_NAME ]
    then
        echo scancel $num
    fi
done
