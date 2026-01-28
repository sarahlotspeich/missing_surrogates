#!/bin/bash
#SBATCH --mail-user=lotspes@wfu.edu  # email address
#SBATCH --mail-type=BEGIN  # Alerts sent when job begins, ends, or aborts
#SBATCH --nodes=1   # Number of nodes required
#SBATCH --ntasks=1   # Number of nodes required
#SBATCH --mem=64G  # Total memory (RAM) required, per node
#SBATCH --time=3-21:30:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]
#SBATCH --array=0-19
#SBATCH --job-name=misspec # job name
#SBATCH --partition=large

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load apps/r/4.3.3 
Rscript sett4_misspec_ipw_smle.R $SLURM_ARRAY_TASK_ID
