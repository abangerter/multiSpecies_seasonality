#! /bin/bash

#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=96:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab

# define file name & argument from array
echo "load in array text"
argFile="/scratch/ab5dr/multiSpecies/scripts/MSlist_noMerge.txt"
#argFile="/scratch/ab5dr/multiSpecies/scripts/MSlist_merge.txt"

echo "define line of array to use as argument"
lineNo=$SLURM_ARRAY_TASK_ID

echo "pick out the sample file ID"
poolID=$(sed -n ${lineNo}p ${argFile})

echo "begin mapping"
./map_reads_merge.sh ${poolID}
