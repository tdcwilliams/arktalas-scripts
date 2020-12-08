#!/bin/bash -x
## Project:
#SBATCH --account=nn2993k
## Job name:
#SBATCH --job-name=ns81-6kPa-h3-tp15d
## Output file
#SBATCH --output=logs/slurm.%j.log   # Stdout & stderr
## Wall time limit:
##SBATCH --time=0-00:20:00
#SBATCH --time=0-12:00:00
## Number of nodes:
#SBATCH --nodes=1
## Number of tasks (total)
##SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
## Set OMP_NUM_THREADS
#SBATCH --cpus-per-task=1

# Queue Option
#  devel: <= 30mins
#  short: <= 120mins
#  normal: > 120mins
#  bigmem?
#SBATCH --qos=preproc
##SBATCH --qos=devel

## Email info
#SBATCH --mail-type=ALL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=timothy.williams@nersc.no # email to the user

function _read_duration {
    tmp=(`grep "duration=" inputs/nextsim.cfg`)
    for tmp_ in ${tmp[@]}
    do
        duration=${tmp_#*=}
        _opt=${tmp_%=*}
        [[ "$_opt" == "duration" ]] && break || unset duration
    done
    [[ -z $duration ]] && { echo "Could not determine duration"; exit 1; }
}

function _check_duration {
    # check running for $split_length days will not take us after $date_final
    # - reduce $tmp_duration if so
    tmp_duration=$split_length
    date_end=`date -d "$START + $tmp_duration days" "+%Y%m%d"`
    while [ $date_end -gt $date_final ]
    do
        (( --tmp_duration ))
        date_end=`date -d "$START + $tmp_duration days" "+%Y%m%d"`
    done
}

function _check_restart_file {
    # check if there is a restart file for start date
    # - if not, either use no restart or use [field,mesh]_final.[bin,dat]
    tmp=`grep "exporter_path=" inputs/nextsim.cfg`
    _restart_path="${tmp#*=}/restart"
    restart_base="${START}T000000Z"
    _restart_file="$_restart_path/field_${restart_base}.bin"
    [[ ! -f $_restart_file ]] && restart_base="final"
}

cd $SLURM_SUBMIT_DIR
source $HOME/pynextsim.sing.src
split_length=20
mpirun_cmd="mpirun tmp/nextsim.exec \
    -mat_mumps_icntl_23 2000 \
    --config-files=inputs/nextsim.cfg"

# start and final day
tmp=`grep "time_init=" inputs/nextsim.cfg`
date_begin=`date -d "${tmp#*=}" "+%Y%m%d"`
_read_duration #read $duration
date_final=`date -d "$date_begin + $duration days" "+%Y%m%d"`
[[ -z $START ]] && { echo "Please enter start date with --export=START=yyyymmdd"; exit 1; }
[[ $START -lt $date_begin ]] && { echo "Bad input start date (<$date_begin)"; exit 1; }
[[ $START -ge $date_final ]] && { echo "Bad input start date (>$date_final)"; exit 1; }
_check_duration # set date_end, tmp_duration

# Special options for first day/restarts
ns_opts=("--simul.duration=$tmp_duration")
[[ $START -eq $date_begin ]] && ns_opts+=(--restart.start_from_restart=false)
_check_restart_file # set restart_base
ns_opts+=(--restart.basename=$restart_base)

# run the model
$mpirun_cmd ${ns_opts[@]} &> "logs/nextsim.${SLURM_JOB_ID}.log" || exit 1

# Only resubmit if we haven't reached $date_final
[[ $date_end -lt $date_final ]] && sbatch --export=START=$date_end $0
exit $?
