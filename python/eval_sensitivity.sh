#! /bin/bash -x
function usage {
    echo $(basename $0) RUNLIST STEP MOORINGS_MASK
    echo RUNLIST is .csv file with the experiment list
    echo STEP is 1 or 2
    echo MOORINGS_MASK is date format of moorings file eg Moorings.nc or Moorings_%Ym%m.nc
    exit 1
}
[[ $# -ne 3 ]] && usage
RUNLIST=$1
STEP=$2
MOORINGS_MASK=$3

[[ $STEP -eq 1 ]] && monthly_eval=1
[[ $STEP -eq 2 ]] && monthly_eval=1

[[ ! -f $RUNLIST ]] && { echo $RUNLIST not found; exit 1; }
runlist=($(cat $RUNLIST))
for run in "${runlist[@]}"
do
    n=${#run}
    run_=${run#*,}
    n_=${#run_}
    edir=$ROOT_DIR/${run:0:$((n-n_-1))}
    [[ "$edir" == "Experiment Directory" ]] && continue
    ./eval_freerun.sh $edir $monthly_eval $MOORINGS_MASK 1 1 1
done
