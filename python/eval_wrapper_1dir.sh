#! /bin/bash -x
function usage {
    echo "$(basename $0) INDIR STEP [MOORINGS_MASK]"
    echo "INDIR is directory with a directory called outputs"
    echo "    and evaluation outputs to go into directories eval-*"
    echo "STEP is 1, 2 or 3"
    echo "    1: do monthly evaluations"
    echo "    2: clean and merge the time-series files/plots"
    echo "    3: collect the monthly maps into pdf files and merge them into 1"
    echo "MOORINGS_MASK is date format of moorings file"
    echo "    eg Moorings.nc (default) or Moorings_%Ym%m.nc"
    exit 1
}
[[ $# -lt 2 ]] && usage
INDIR=$1
STEP=$2
MOORINGS_MASK=${3-"Moorings.nc"}

if [[ $STEP -eq 3 ]]
then
    singularity exec --cleanenv $PYNEXTSIM_SIF \
        ./collect_maps.sh $INDIR
    exit $?
fi

[[ $STEP -eq 1 ]] && monthly_eval=1
[[ $STEP -eq 2 ]] && monthly_eval=0
./eval_freerun.sh $INDIR $monthly_eval $MOORINGS_MASK 1 1 1
