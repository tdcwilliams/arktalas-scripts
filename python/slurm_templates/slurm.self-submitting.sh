#!/bin/bash -x
#SBATCH --account=%%{proj_num}
#SBATCH --job-name=%%{job_name}
#SBATCH --output=logs/slurm.%j.log   # Stdout & stderr
#SBATCH --time=%%{wall_time}
#SBATCH --nodes=%%{num_nodes}
#SBATCH --ntasks=%%{num_tasks}
#SBATCH --cpus-per-task=1

## Email info
#SBATCH --mail-type=ALL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=%%{email}

SPLIT_LENGTH=%%{split_length}
source %%{env_file}

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
    # check running for $SPLIT_LENGTH days will not take us after $date_final
    # - reduce $tmp_duration if so
    tmp_duration=$SPLIT_LENGTH
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
    restart_path="${tmp#*=}/restart"
    restart_base="${START}T000000Z"
    _restart_file="$restart_path/field_${restart_base}.bin"
    [[ ! -f $_restart_file ]] && unset restart_base
}

function _check_start_date {
    # check start date is allowed
    [[ $START -lt $date_begin ]] && \
        { echo "Bad input start date (<$date_begin)"; exit 1; }
    [[ $START -ge $date_final ]] &&  \
        { echo "Bad input start date (>$date_final)"; exit 1; }
}

# init
cd $SLURM_SUBMIT_DIR
ns_opts=()

# get $duration
_read_duration
tmp=`grep "time_init=" inputs/nextsim.cfg`
date_begin=`date -d "${tmp#*=}" "+%Y%m%d"`
date_final=`date -d "$date_begin + $duration days" "+%Y%m%d"`

# get $START and check it
START=${START-"$date_begin"}
_check_start_date

# set date_end, tmp_duration
_check_duration
ns_opts+=("--simul.duration=$tmp_duration")

# set restart_base and restart_path
# - default is set in nextsim.cfg
_check_restart_file
[[ ! -z $restart_base ]] && ns_opts+=(--restart.basename=$restart_base)
if [ $START -gt $date_begin ]
then
    # Set extra options if not first day
    ns_opts+=(--simul.spinup_duration=0)
    ns_opts+=(--restart.start_from_restart=true)
    ns_opts+=(--restart.input_path=$restart_path)
    # assume we start from final if not on first day and restart base hasn't been set already
    [[ -z $restart_base ]] && ns_opts+=(--restart.basename=final)
fi

# run the model
mpirun -np $SLURM_NTASKS tmp/nextsim.exec \
    --config-files=inputs/nextsim.cfg ${ns_opts[@]} \
    &> "logs/nextsim.${SLURM_JOB_ID}.log" || exit 1

# Exit if we have reached $date_final
[[ $date_end -eq $date_final ]] && exit 0

# Otherwise resubmit and exit
sbatch --export=START=$date_end $0
exit $?
