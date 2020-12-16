#! /bin/bash -x

if [ $# -lt 3 ]
then
    echo "Usage"
    echo "eval_freerun.sh INPUT_DIR MONTHLY_EVAL MOORINGS_MASK DO_CS2SMOS [DO_OSISAF] [DO_DRIFT] [DO_PLOTS] [DO_SMOS] [DO_AMSR2]"
    exit 0
fi
ME=`readlink -f $0`
HERE=`dirname $ME`

do_copy=0 # copy moorings if running
FCDIR=$1
MONTHLY_EVAL=$2
MOORINGS_MASK=${3-"Moorings_%Ym%m.nc"}
DO_CS2SMOS=${4}
DO_OSISAF=${5-"0"}
DO_DRIFT=${6-"0"}
DO_PLOTS=${7-"0"}
DO_SMOS=${8-"0"}
DO_AMSR2=${9-"0"}
numproc_sem=12

function run
{
    echo $1    

    # if not on fram, launch normally
    #[[ "${HOSTNAME:10:14}" != "fram.sigma2.no" ]] && { $1; return; }
    [[ "${HOSTNAME:10:14}" != "fram.sigma2.no" ]] && { sem -j $numproc_sem $1 & return; }

    # on fram, launch with sbatch
    mkdir -p logs
    args=($1)
    sbatch_opts=()
    sbatch_opts+=("--job-name=AN${args[0]}")
    sbatch_opts+=("--time=02:30:00")
    sbatch \
        ${sbatch_opts[@]} \
        --export=COMMAND="singularity exec --cleanenv $PYNEXTSIM_SIF $1" \
        $HERE/slurm_batch_job.sh
}
function take_monthly_averages
{
    # take monthly averages and save hi-res pics
    # - can skip months with $SKIP
    outdir=$1
    [[ -z $outdir ]] && { echo "no output dir provided"; exit 1; }
    for mon in `seq 0 $N_MONTHS`
    do
        date0=`date -d "${MONTH0}01 + $mon months" "+%Y%m%d"`
        [[ "${SKIP[@]}" =~ `date -d "$date0" "+%m"` ]] && continue
    
        date1=`date -d "$date0 + 1 month" "+%Y%m%d"`
        date1=`date -d "$date1 - 1 day" "+%Y%m%d"`
        run "$CMD -d0 $date0 -d1 $date1 -o $outdir"
    done
}
function do_standard_run
{
    # do eval but don't plot maps
    outdir=$1
    [[ -z $outdir ]] && { echo "no output dir provided"; exit 1; }
    rm -f $outdir/*.txt $outdir/*.png #clean old time series files and figures
    for lim in "${LIMS[@]}"
    do
        run "$CMD $lim -o $outdir -nm"
    done
}

#copy Moorings to avoid causing a run to crash
fcdir=$FCDIR/outputs

# Inputs that are independant of the source
inputs="$fcdir -np -g medium_arctic_10km.msh"
if [ $MONTHLY_EVAL -eq 1 ]
then
    MONTH0=201211
    N_MONTHS=12
fi

if [ $DO_CS2SMOS -eq 1 ]
then
   smoothing="-sig 1"
   CMD="evaluate_forecast.py $inputs $smoothing -s Cs2SmosThick -mm $MOORINGS_MASK"
   LIMS=("-d1 20130430" "-d0 20131015")
   SKIP=(05 06 07 08 09) #skip averages of "summer" months
   odir="$FCDIR/eval-cs2smos"
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" \
       || do_standard_run "$odir"
fi

if [ $DO_OSISAF -eq 1 ]
then
   smoothing="-sig 2"
   CMD="evaluate_forecast.py $inputs $smoothing -s OsisafConc -mm $MOORINGS_MASK"
   LIMS=("")
   SKIP=() #don't skip any monthly averages
   odir="$FCDIR/eval-osisaf-conc"
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" || \
       do_standard_run "$odir"

   # blend AMSR2
   CMD+=" -b"
   odir="$FCDIR/eval-osisaf-amsr2-conc"
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" || \
       do_standard_run "$odir"

   # extent
   CMD="evaluate_forecast.py $inputs $smoothing -s OsisafConc -ee -mm $MOORINGS_MASK"
   odir=$FCDIR/eval-osisaf-extent
   [[ $MONTHLY_EVAL -eq 1 ]] && CMD+=" -nm" # don't need maps
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" \
       || do_standard_run "$odir"

   # blend AMSR2
   CMD+=" -b"
   odir="$FCDIR/eval-osisaf-amsr2-extent"
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" || \
       do_standard_run "$odir"
fi

if [ $DO_DRIFT -eq 1 ]
then
   inputs_drift="$fcdir -mu 2.5"
   #inputs_drift+=" -f"
   #inputs_drift="$FCDIR -np -g medium_arctic_10km.msh"
   CMD="evaluate_drift_forecast.py $inputs_drift -mm $MOORINGS_MASK"
   LIMS=("-d1 20130430") # finish 1st eval round after 2019-04
   LIMS+=("-d0 20131001") # start 2nd eval round at 2019-10
   SKIP=(05 06 07 08 09) #skip averages of "summer" months
   odir="$FCDIR/eval-osisaf-drift"
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" \
       || do_standard_run "$odir"

   inputs_drift="$fcdir -mu 20"
   #inputs_drift+=" -f"
   #inputs_drift="$FCDIR -np -g medium_arctic_10km.msh"
   CMD="evaluate_drift_forecast.py $inputs_drift -mm $MOORINGS_MASK"
   LIMS=("")
   SKIP=() #include "summer" months
   odir="$FCDIR/eval-osisaf-drift-mu10kpd"
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$odir" \
       || do_standard_run "$odir"
fi

if [ $DO_PLOTS -eq 1 ]
then
    CMD="plot_nextsim_output.py $FCDIR plot.cfg -o $FCDIR/figs"
    LIMS=("")
    do_standard_run
fi

if [ $DO_SMOS -eq 1 ]
then
   smoothing="-sf 0.5"
   #smoothing="-sig 1"
   CMD="evaluate_forecast.py $inputs -nb 1 $smoothing -s SmosThick -mm $MOORINGS_MASK"
   LIMS=("-d1 20190430" "-d0 20191015")
   SKIP=(05 06 07 08 09) #skip averages of "summer" months
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$FCDIR/eval-smos" \
       || do_standard_run "$FCDIR/eval-smos"
fi

if [ $DO_AMSR2 -eq 1 ]
then
   CMD="evaluate_forecast.py $inputs -s Amsr2Conc -mm $MOORINGS_MASK"
   LIMS=("")
   SKIP=() #don't skip any monthly averages
   [[ $MONTHLY_EVAL -eq 1 ]] \
       && take_monthly_averages "$FCDIR/eval-amsr2" \
       || do_standard_run "$FCDIR/eval-amsr2"
fi
