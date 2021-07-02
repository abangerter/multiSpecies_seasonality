#! /bin/bash

#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab


### set paths
	VarDir="/home/ab5dr/poolSeqPipeline/varScan"
	outDir="/scratch/ab5dr/multiSpecies/deNovoCalls"
	species=${1}

### Do VarScan for indels
	echo "run VarScan"
	java -Xmx40g -jar ${VarDir}/VarScan.v2.3.9.jar mpileup2indel ${outDir}/MS_${species}.forVarScan.mpileup --output-vcf \
	> ${outDir}/MS_${species}.VarScan.indel.vcf
