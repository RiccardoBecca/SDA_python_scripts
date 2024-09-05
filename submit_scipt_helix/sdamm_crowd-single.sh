#!/bin/bash
#SBATCH -e tryp-BEN-60C-%j.err
#SBATCH -o tryp-BEN-60C-%j.out
#SBATCH --partition=single
#SBATCH --time=50:00:00
##SBATCH -ntasks-per-node=1
##SBATCH --ntasks=8
##SBATCH --exclusive
#SBATCH --cpus-per-task=8
##SBATCH --mem-per-cpu=2G
#SBATCH --array=1-1000

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SDA_FLEX=/home/hits/hits_hits/hits_beccarro/sda-git-workstation

## try for sda_association

#START=$((${SLURM_ARRAY_TASK_ID}*${SLURM_NTASKS}+1))
#END=$((${SLURM_ARRAY_TASK_ID}*${SLURM_NTASKS}+${SLURM_NTASKS}))

#for (( i=${START}; i<=${END}; i++ ))
#do
##check if file exist
$SDA_FLEX/bin/sda_flex sdamm_crowd_${SLURM_ARRAY_TASK_ID}.in > sdamm_${SLURM_ARRAY_TASK_ID}.out
#done


echo "Time in seconds is: "$SECONDS
