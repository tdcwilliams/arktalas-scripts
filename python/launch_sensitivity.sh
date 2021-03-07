#! /bin/bash -x
function usage {
    echo Usage:
    echo "$(basename $0) ROOT_DIR BATCH_NAME CONFIG_FILE NS_CONFIG_FILE EXECUTABLE [LAUNCH]"
    exit 1
}
[[ $# -lt 5 ]] && usage
ROOT_DIR=$1
BATCH_NAME=$2
CONFIG_FILE=$3
NS_CONFIG_FILE=$4
EXECUTABLE=$5
LAUNCH=${6-1}

# setup
singularity exec --cleanenv $PYNEXTSIM_SIF \
    ./setup_sensitivity.py $ROOT_DIR $BATCH_NAME $CONFIG_FILE $NS_CONFIG_FILE $EXECUTABLE \
    || exit 1

# launch if we want to
runlist=$ROOT_DIR/$BATCH_NAME.csv
msg="Run summary in $runlist"
[[ LAUNCH -eq 0 ]] && { echo ""; echo $msg; exit 0; }

[[ ! -f $runlist ]] && { echo $runlist not found; exit 1; }
runlist=($(cat $runlist))
for run in "${runlist[@]}"
do
    n=${#run}
    run_=${run#*,}
    n_=${#run_}
    edir=$ROOT_DIR/${run:0:$((n-n_-1))}
    [[ "$edir" == "Experiment Directory" ]] && continue
    cd $edir && sbatch inputs/slurm.sh
done

echo ""
echo $msg
