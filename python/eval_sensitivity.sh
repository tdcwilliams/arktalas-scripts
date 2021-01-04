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

# loop over all the dirs in the .csv file
# and put them through eval_wrapper_1dir.sh
[[ ! -f $RUNLIST ]] && { echo $RUNLIST not found; exit 1; }
runlist=($(cat $RUNLIST |grep -v Experiment))
rootdir=$(dirname $RUNLIST)
for run in "${runlist[@]}"
do
    n=${#run}
    run_=${run#*,}
    n_=${#run_}
    edir=$rootdir/${run:0:$((n-n_-1))}
    ./eval_wrapper_1dir.sh $edir $STEP $MOORINGS_MASK
done

[[ $STEP != 3 ]] && exit 0
# compare scalar and drift metrics
batch_name=$(basename $RUNLIST)
batch_name=${batch_name%.csv}
batch_name=${batch_name#sens_}
odir=figs/comp_runs/comp_runs.$batch_name
for t in scalars drift
do
    s=comp_runs_$t
    singularity exec --cleanenv $PYNEXTSIM_SIF \
        ./${s}.py config-files/${s}.${batch_name}.cfg -o $odir
done
