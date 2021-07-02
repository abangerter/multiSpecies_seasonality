#! /bin/bash


### modified by Alyssa Bangerter from a script written by Martin Kapun on 10/09/15
### modified on 03/22/2021 -- updated to new ref genomes for Dsim and Zind 
### modified last on 4/19/2021 -- updated to a newer Dsuz ref genome


### set incoming arguments
   echo "${1}"
   name=${1} # name of the output-file, e.g. Barcelona-spring
   flowcell=${2}
   read1=${flowcell}_s2_1_${name}.fastq.gz # full path to Read1
   read2=${flowcell}_s2_2_${name}.fastq.gz # full path to Read2


### set other variable paths
   echo "setting paths & variables"
   samtoolsDir="/home/ab5dr/poolSeqPipeline/samtools/samtools-1.3"
   cutDir="/home/ab5dr/poolSeqPipeline/cutadapt/cutadapt-1.8.3/bin"
   inDir="/project/berglandlab/alyssa/novaSeqRaw"
   fastQCdir="/home/ab5dr/poolSeqPipeline/fastQC/FastQC.app/Contents/Resources/Java"
   outDir="/scratch/ab5dr/multiSpecies"
   finalDir="/scratch/ab5dr/multiSpecies/realignBam"
   picardDir="/home/ab5dr/poolSeqPipeline/picard/picard-tools-1.109"
   gatkDir="/home/ab5dr/poolSeqPipeline/gatk"


### select the reference genome 
   species=($( echo ${name} | cut -f 4 -d "_" ))
   if [ ${species} = "Dmel" ]; then
      ref="/scratch/ab5dr/pooledSeq2019/poolSeqPipe/misc/holo_dmel_6.12.fa"
   elif [ ${species} = "Dsim" ]; then
      ref="/scratch/ab5dr/multiSpecies/ref/GCF_016746395.1_Prin_Dsim_3.0_genomic.fna"
   elif [ ${species} = "Dhyd" ]; then
      ref="/scratch/ab5dr/multiSpecies/ref/GCF_003285905.1_DhydRS2_genomic.fna"
   elif [ ${species} = "Zind" ]; then  
      ref="/scratch/ab5dr/multiSpecies/ref/z_indianus_16GNV01_v02.fasta"
   elif [ ${species} = "Dsuz" ]; then
      ref="/scratch/ab5dr/multiSpecies/ref/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.fna"
   else
      echo "reference determination error"
   fi


### for Rivanna, load modules for needed programs
   echo "loading modules"
   module load gcc/7.1.0
   module load bwa/0.7.15


### run fastQC 
   echo "fastQC"
   perl ${fastQCdir}/fastqc ${inDir}/${read1} ${inDir}/${read2} -o ${outDir}/fastQC/raw


### trim out adapter sequences 
   echo "trimming out adapter sequences"
   ${cutDir}/cutadapt \
      -q 18 \
      --minimum-length 75 \
      -o ${outDir}/trimmedFQ/${name}_1.fq.gz \
      -p ${outDir}/trimmedFQ/${name}_2.fq.gz \
      -b  GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
      -b  GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
      -O 15 \
      -n 3 \
      ${inDir}/${read1} \
      ${inDir}/${read2}
#-b ACACTCTTTCCCTACACGACGCTCTTCCGATC
#-B CAAGCAGAAGACGGCATACGAGAT
#-b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#-b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#-b GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#-b CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
#-b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
#-b GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG


### fastQC on the trimmed files
   echo "fastQC on trimmed files"
   perl ${fastQCdir}/fastqc ${outDir}/trimmedFQ/${name}_1.fq.gz ${outDir}/trimmedFQ/${name}_2.fq.gz -o ${outDir}/fastQC/trim


### map with BWA
   echo "mapping"
   bwa mem -M -t $SLURM_CPUS_PER_TASK \
      -R "@RG\tID:${name}\tLB:${name}\tPL:illumina\tSM:sample\tPU:${name}" \
      ${ref} \
      ${outDir}/trimmedFQ/${name}_1.fq.gz \
      ${outDir}/trimmedFQ/${name}_2.fq.gz | \
   ${samtoolsDir}/samtools view -Sbh -q 20 -F 0x100 > ${outDir}/mappedBam/${name}.bam


### sort mapped output
   echo "sorting"
   java \
   -Xmx20g \
   -Dsnappy.disable=true \
   -Djava.io.tmpdir=${outDir}/tmp \
   -jar ${picardDir}/SortSam.jar \
   I=${outDir}/mappedBam/${name}.bam \
   O=${outDir}/sortBam/${name}.sort.bam \
   SO=coordinate \
   VALIDATION_STRINGENCY=SILENT


