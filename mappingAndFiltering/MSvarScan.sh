#! /bin/bash

#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=96:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab


### set paths
	samtoolsDir="/home/ab5dr/poolSeqPipeline/samtools/samtools-1.3"
	VarDir="/home/ab5dr/poolSeqPipeline/varScan"
	outDir="/scratch/ab5dr/multiSpecies/deNovoCalls"


### select the reference genome and bamlist
   species=${1}
   if [ ${species} = "Dmel" ]; then
      ref="/scratch/ab5dr/pooledSeq2019/poolSeqPipe/misc/holo_dmel_6.12.fa"
      bamlist="/scratch/ab5dr/multiSpecies/scripts/MS_Dmel.bamlist"
   elif [ ${species} = "Dsim" ]; then
      ref="/scratch/ab5dr/multiSpecies/ref/GCF_016746395.1_Prin_Dsim_3.0_genomic.fna"
      bamlist="/scratch/ab5dr/multiSpecies/scripts/MS_Dsim.bamlist"
   elif [ ${species} = "Dhyd" ]; then
      ref="/scratch/ab5dr/multiSpecies/ref/GCF_003285905.1_DhydRS2_genomic.fna"
      bamlist="/scratch/ab5dr/multiSpecies/scripts/MS_Dhyd.bamlist"
   elif [ ${species} = "Zind" ]; then  
      ref="/scratch/ab5dr/multiSpecies/ref/z_indianus_16GNV01_v02.fasta"
      bamlist="/scratch/ab5dr/multiSpecies/scripts/MS_Zind.bamlist"
   elif [ ${species} = "Dsuz" ]; then
      ref="/scratch/ab5dr/multiSpecies/ref/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.fna"
      bamlist="/scratch/ab5dr/multiSpecies/scripts/MS_Dsuz.bamlist"
   else
      echo "error"
   fi



### Do the mpileup
	echo "run mpileup"
	${samtoolsDir}/samtools mpileup \
	-q 20 -Q 20 -B \
	-b ${bamlist} \
	-f ${ref} \
	> ${outDir}/MS_${species}.forVarScan.mpileup


### Do VarScan
	echo "run VarScan"
	java -Xmx40g -jar ${VarDir}/VarScan.v2.3.9.jar mpileup2snp ${outDir}/MS_${species}.forVarScan.mpileup --min-var-freq 0.005 --min-avg-qual 20 --output-vcf --variants \
	> ${outDir}/MS_${species}.VarScan.vcf


