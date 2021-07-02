#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
##SBATCH --mem=64000
#SBATCH --time=96:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

# define file name & argument from array
echo "load in array text"
argFile="/scratch/ab5dr/multiSpecies/scripts/MSlist.txt"

echo "define line of array to use as argument"
lineNo=$SLURM_ARRAY_TASK_ID
export SLURM_CPUS_PER_TASK

echo "pick out the sample file ID"
args=$(sed -n ${lineNo}p ${argFile})

echo "begin mapping"
./map_reads.sh ${args}
