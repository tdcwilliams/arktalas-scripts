#! /bin/bash -x
#SBATCH --account=nn2993k
#SBATCH --job-name=JOBNAME
#SBATCH --time 30:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --output=logs/slurm.%j.out

function usage {
    echo "sbatch --export=COMMAND=\"command\",SRC_FILE=\"path-to-source-file\" `basename $0`"
}
[[ -z $COMMAND ]] && { usage; exit 1; }
SRC_FILE=${SRC_FILE-"$HOME/pynextsim.sing.src"}
source $SRC_FILE
cd $SLURM_SUBMIT_DIR
$COMMAND
