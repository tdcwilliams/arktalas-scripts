#!/bin/bash -x
#SBATCH --account=%%{PROJ_NUM}
#SBATCH --job-name=%%{JOB_NAME}
#SBATCH --output=logs/slurm.%j.log
#SBATCH --time=%%{WALL_TIME}
#SBATCH --nodes=%%{NUM_NODES}
#SBATCH --ntasks=%%{NUM_TASKS}
#SBATCH --cpus-per-task=1
#SBATCH --qos=%%{QOS}
##SBATCH --partition=bigmem
##SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%%{EMAIL}

cd $SLURM_SUBMIT_DIR
source %%{ENV_FILE}
singularity exec --cleanenv $NEXTSIM_SIF \
    mpirun tmp/nextsim.exec --config-files=inputs/nextsim.cfg \
    &> "logs/nextsim.${SLURM_JOB_ID}.log" || exit 1
