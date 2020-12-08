#!/bin/bash -x
#SBATCH --account=nn2993k
#SBATCH --job-name=breakup-00-wrf10km
#SBATCH --output=logs/slurm.%j.log
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
##SBATCH --partition=bigmem
##SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothy.williams@nersc.no

cd $SLURM_SUBMIT_DIR
source $HOME/pynextsim.sing.src
singularity exec --cleanenv $NEXTSIM_SIF \
    mpirun tmp/nextsim.exec --config-files=inputs/nextsim.cfg \
    &> "logs/nextsim.${SLURM_JOB_ID}.log" || exit 1
