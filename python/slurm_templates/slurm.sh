#!/bin/bash -x
#SBATCH --account=%%{proj_num}
#SBATCH --job-name=%%{job_name}
#SBATCH --output=logs/slurm.%j.log
#SBATCH --time=%%{wall_time}
#SBATCH --nodes=%%{num_nodes}
#SBATCH --ntasks=%%{num_tasks}
#SBATCH --cpus-per-task=1
#%%{qos_line}
#%%{bigmem_line}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%%{email}

cd $SLURM_SUBMIT_DIR
source %%{env_file}
mpirun tmp/nextsim.exec --config-files=inputs/nextsim.cfg \
    &> "logs/nextsim.${SLURM_JOB_ID}.log" || exit 1
