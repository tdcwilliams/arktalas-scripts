#! /bin/bash -x
function usage {
    echo "$(basename $0) RUNLIST STEP [MOORINGS_MASK]"
    echo "RUNLIST is .csv file with the experiment list"
    echo "STEP is 1, 2 or 3"
    echo "    1: do monthly evaluations"
    echo "    2: clean and merge the time-series files/plots"
    echo "    3: collect the monthly maps into pdf files and merge them into 1"
    echo "MOORINGS_MASK is date format of moorings file"
    echo "    eg Moorings.nc (default) or Moorings_%Ym%m.nc"
    exit 1
}
[[ $# -lt 2 ]] && usage
RUNLIST=$1
STEP=$2
MOORINGS_MASK=${3-Moorings.nc}

[[ $STEP -eq 1 ]] && monthly_eval=1
[[ $STEP -eq 2 ]] && monthly_eval=1
[[ $STEP -eq 3 ]] && \
    collect="singularity exec --cleanenv $PYNEXTSIM_SIF ./collect_maps.sh"

[[ ! -f $RUNLIST ]] && { echo $RUNLIST not found; exit 1; }
runlist=($(cat $RUNLIST))
for run in "${runlist[@]}"
do
    n=${#run}
    run_=${run#*,}
    n_=${#run_}
    edir=$ROOT_DIR/${run:0:$((n-n_-1))}
    [[ "$edir" == "Experiment Directory" ]] && continue
    [[ $STEP -eq 3 ]] && $collect $edir && continue
    ./eval_freerun.sh $edir $monthly_eval $MOORINGS_MASK 1 1 1
done
