#! /bin/bash


### modified by Alyssa Bangerter from a script written by Martin Kapun on 01/10/15 to work on multiple Drosophilids 
### modified on 03/22/2021 -- edited to include new Dsim and Zind reference genomes 
### modified last on 04/29/2021 -- edited to include newer Dsuz reference genome
### find python scripts in separate directory on this github


### set incoming arguments
	echo "${1}"
	name=${1} # name of the sample


### set other variable paths
	echo "setting paths & variables"
	samtoolsDir="/home/ab5dr/poolSeqPipeline/samtools/samtools-1.3"
	bamFQdir="/home/ab5dr/poolSeqPipeline/bam2fastq/bam2fastq-1.1.0"
	inDir="/scratch/ab5dr/multiSpecies/realignBam"
	interDir="/scratch/ab5dr/multiSpecies"
	outDir="/scratch/ab5dr/multiSpecies/poolSeqDecontam"


### for Rivanna, load modules for needed programs
	echo "loading modules"
	module load gcc/7.1.0
	module load bwa/0.7.15
	module load openmpi/3.1.4
	module load python/2.7.16


### turn mapped and realigned bam files into fastq files
	echo "turning bams into fastq"
	${bamFQdir}/bam2fastq -s -o ${interDir}/bam2fastqOut/${name}.dedup.indel# ${inDir}/${name}.dedup.indel.bam
	gzip ${interDir}/bam2fastqOut/${name}.dedup.indel*


### make ref genomes that combine the different genomes - run once, comment out afterwards
	# mel as the contaminant
		#sed 's/>/>mel_/g' /scratch/ab5dr/pooledSeq2019/poolSeqPipe/misc/holo_dmel_6.12.fa | gzip -c > /scratch/ab5dr/multiSpecies/ref/mel_genome_prefix.fa.gz
		#zcat /scratch/ab5dr/multiSpecies/ref/mel_genome_prefix.fa.gz | cat /scratch/ab5dr/multiSpecies/ref/GCF_016746395.1_Prin_Dsim_3.0_genomic.fna - | gzip -c > /scratch/ab5dr/multiSpecies/ref/melContam_sim_combo.fa.gz
		#bwa index -a bwtsw /scratch/ab5dr/multiSpecies/ref/melContam_sim_combo.fa.gz
	# sim as contaminant
		#sed 's/>/>sim_/g' /scratch/ab5dr/multiSpecies/ref/GCF_016746395.1_Prin_Dsim_3.0_genomic.fna | gzip -c > /scratch/ab5dr/multiSpecies/ref/sim_genome_prefix.fa.gz
		#zcat /scratch/ab5dr/multiSpecies/ref/sim_genome_prefix.fa.gz | cat /scratch/ab5dr/pooledSeq2019/poolSeqPipe/misc/holo_dmel_6.12.fa - | gzip -c > /scratch/ab5dr/multiSpecies/ref/simContam_mel_combo.fa.gz
		#bwa index -a bwtsw /scratch/ab5dr/multiSpecies/ref/simContam_mel_combo.fa.gz
	# mel/sim as contam 
		#cat /scratch/ab5dr/pooledSeq2019/poolSeqPipe/misc/holo_dmel_6.12.fa | cat /scratch/ab5dr/multiSpecies/ref/GCF_016746395.1_Prin_Dsim_3.0_genomic.fna > /scratch/ab5dr/multiSpecies/ref/mel_sim_combo.fa
		#sed 's/>/>MelSim_/g' /scratch/ab5dr/multiSpecies/ref/mel_sim_combo.fa | gzip -c > /scratch/ab5dr/multiSpecies/ref/combined_melsim_prefix.fa.gz
	# D. hyd.
		#zcat /scratch/ab5dr/multiSpecies/ref/combined_melsim_prefix.fa.gz | cat /scratch/ab5dr/multiSpecies/ref/GCF_003285905.1_DhydRS2_genomic.fna - | gzip -c > /scratch/ab5dr/multiSpecies/ref/hyd_mel_sim_combo.fa.gz
		#bwa index -a bwtsw /scratch/ab5dr/multiSpecies/ref/hyd_mel_sim_combo.fa.gz
	# Z. ind.
		#zcat /scratch/ab5dr/multiSpecies/ref/combined_melsim_prefix.fa.gz | cat /scratch/ab5dr/multiSpecies/ref/z_indianus_16GNV01_v02.fasta - | gzip -c > /scratch/ab5dr/multiSpecies/ref/ind_mel_sim_combo.fa.gz
		#bwa index -a bwtsw /scratch/ab5dr/multiSpecies/ref/ind_mel_sim_combo.fa.gz
	# D. suz.
		#zcat /scratch/ab5dr/multiSpecies/ref/combined_melsim_prefix.fa.gz | cat /scratch/ab5dr/multiSpecies/ref/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.fna - | gzip -c > /scratch/ab5dr/multiSpecies/ref/suz_mel_sim_combo.fa.gz
		#bwa index -a bwtsw /scratch/ab5dr/multiSpecies/ref/suz_mel_sim_combo.fa.gz


### species-based run  
   species=($( echo ${name} | cut -f 4 -d "_" ))
   if [ ${species} = "Dmel" ]; then
   		holoref="/scratch/ab5dr/multiSpecies/ref/simContam_mel_combo.fa.gz"
   		### competitive mapping 
			echo "competitive mapping"
			bwa mem -Mt $SLURM_CPUS_PER_TASK ${holoref} ${interDir}/bam2fastqOut/${name}.dedup.indel\_1.gz ${interDir}/bam2fastqOut/${name}.dedup.indel\_2.gz > ${interDir}/compMap/${name}.dedup.indel.compMap.sam
		### de-convolute  
			echo "de-convolute"
			export PATH=$PATH:/home/ab5dr/poolSeqPipeline/samtools/samtools-0.1.19
			python ${interDir}/scripts/fix_bam.py \
			--contaminated ${inDir}/${name}.dedup.indel.bam \
			--prefix sim_ \
			--detect ${interDir}/compMap/${name}.dedup.indel.compMap.sam \
			--output ${outDir}/${name}.dedup.indel.compMap
   elif [ ${species} = "Dsim" ]; then
   		holoref="/scratch/ab5dr/multiSpecies/ref/melContam_sim_combo.fa.gz"
   		### competitive mapping 
			echo "competitive mapping"
			bwa mem -Mt $SLURM_CPUS_PER_TASK ${holoref} ${interDir}/bam2fastqOut/${name}.dedup.indel\_1.gz ${interDir}/bam2fastqOut/${name}.dedup.indel\_2.gz > ${interDir}/compMap/${name}.dedup.indel.compMap.sam
		### de-convolute  
			echo "de-convolute"
			export PATH=$PATH:/home/ab5dr/poolSeqPipeline/samtools/samtools-0.1.19
			python ${interDir}/scripts/fix_bam_Dsim.py \
			--contaminated ${inDir}/${name}.dedup.indel.bam \
			--prefix mel_ \
			--detect ${interDir}/compMap/${name}.dedup.indel.compMap.sam \
			--output ${outDir}/${name}.dedup.indel.compMap
   elif [ ${species} = "Dhyd" ]; then
  		holoref="/scratch/ab5dr/multiSpecies/ref/hyd_mel_sim_combo.fa.gz"
   		### competitive mapping 
			echo "competitive mapping"
			bwa mem -Mt $SLURM_CPUS_PER_TASK ${holoref} ${interDir}/bam2fastqOut/${name}.dedup.indel\_1.gz ${interDir}/bam2fastqOut/${name}.dedup.indel\_2.gz > ${interDir}/compMap/${name}.dedup.indel.compMap.sam
		### de-convolute  
			echo "de-convolute"
			export PATH=$PATH:/home/ab5dr/poolSeqPipeline/samtools/samtools-0.1.19
			python ${interDir}/scripts/fix_bam_Dhyd.py \
			--contaminated ${inDir}/${name}.dedup.indel.bam \
			--prefix MelSim_ \
			--detect ${interDir}/compMap/${name}.dedup.indel.compMap.sam \
			--output ${outDir}/${name}.dedup.indel.compMap
   elif [ ${species} = "Zind" ]; then
  		holoref="/scratch/ab5dr/multiSpecies/ref/ind_mel_sim_combo.fa.gz"
   		### competitive mapping 
			echo "competitive mapping"
			bwa mem -Mt $SLURM_CPUS_PER_TASK ${holoref} ${interDir}/bam2fastqOut/${name}.dedup.indel\_1.gz ${interDir}/bam2fastqOut/${name}.dedup.indel\_2.gz > ${interDir}/compMap/${name}.dedup.indel.compMap.sam
		### de-convolute  
			echo "de-convolute"
			export PATH=$PATH:/home/ab5dr/poolSeqPipeline/samtools/samtools-0.1.19
			python ${interDir}/scripts/fix_bam_Zind.py \
			--contaminated ${inDir}/${name}.dedup.indel.bam \
			--prefix MelSim_ \
			--detect ${interDir}/compMap/${name}.dedup.indel.compMap.sam \
			--output ${outDir}/${name}.dedup.indel.compMap
   elif [ ${species} = "Dsuz" ]; then
  		holoref="/scratch/ab5dr/multiSpecies/ref/suz_mel_sim_combo.fa.gz"
   		### competitive mapping 
			echo "competitive mapping"
			bwa mem -Mt $SLURM_CPUS_PER_TASK ${holoref} ${interDir}/bam2fastqOut/${name}.dedup.indel\_1.gz ${interDir}/bam2fastqOut/${name}.dedup.indel\_2.gz > ${interDir}/compMap/${name}.dedup.indel.compMap.sam
		### de-convolute  
			echo "de-convolute"
			export PATH=$PATH:/home/ab5dr/poolSeqPipeline/samtools/samtools-0.1.19
			python ${interDir}/scripts/fix_bam_Dsuz.py \
			--contaminated ${inDir}/${name}.dedup.indel.bam \
			--prefix MelSim_ \
			--detect ${interDir}/compMap/${name}.dedup.indel.compMap.sam \
			--output ${outDir}/${name}.dedup.indel.compMap
   else
      echo "species determination error"
   fi


