#! /bin/bash



### modified by Alyssa Bangerter from a script written by Martin Kapun on 10/09/15
### modified on 03/22/2021 -- updated to new ref genomes for Dsim and Zind
### modified last on 04/19/2021 -- updated to a newer ref genome for Dsuz
### see notes below on how to modify the script depending on if you run with the MSlist_merge.txt versus the MSlist_noMerge.txt as input


### set incoming arguments
   echo "${1}"
   name=${1} # name of the output-file, e.g. Barcelona-spring


### set other variable paths
   echo "setting paths & variables"
   samtoolsDir="/home/ab5dr/poolSeqPipeline/samtools/samtools-1.3"
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


### edit v2 read groups 
   ## comment these lines out when running the no merge list
   echo "editing read groups"
   java \
   -Xmx20g \
   -Dsnappy.disable=true \
   -Djava.io.tmpdir=${outDir}/tmp \
   -jar ${picardDir}/AddOrReplaceReadGroups.jar \
   I=${outDir}/sortBam/${name}_v2.sort.bam \
   O=${outDir}/rgSort/${name}_v2.sort.bam \
   RGID=${name}_v2 \
   RGLB=${name} \
   RGPL=illumina \
   RGPU=${name}_v2 \
   RGSM=sample \
   CREATE_INDEX=true \
   VALIDATION_STRINGENCY=SILENT


### merge lanes 1 and 2 data 
   ## comment these lines out when running the no merge list 
   echo "merging lanes"
   mv ${outDir}/sortBam/${name}.sort.bam ${outDir}/rgSort/${name}.sort.bam
   ls ${outDir}/rgSort/${name}*.sort.bam > ${outDir}/listDir/${name}.txt
   ${samtoolsDir}/samtools merge -f -b ${outDir}/listDir/${name}.txt ${outDir}/sortBam/${name}.sort.merge.bam
   ${samtoolsDir}/samtools index ${outDir}/sortBam/${name}.sort.merge.bam


### remove duplicates
   ## swap out input file with other option if running with the no merge list 
   echo "removing duplicates"
   java \
   -Xmx20g \
   -Dsnappy.disable=true \
   -Djava.io.tmpdir=${outDir}/tmp \
   -jar ${picardDir}/MarkDuplicates.jar \
   REMOVE_DUPLICATES=true \
   I=${outDir}/sortBam/${name}.sort.merge.bam \
   O=${outDir}/dedupBam/${name}.dedup.bam \
   M=${outDir}/dedupBam/dedupStats/${name}.txt \
   CREATE_INDEX=true
   VALIDATION_STRINGENCY=SILENT
#   I=${outDir}/sortBam/${name}.sort.bam \


### re-align around indels using GATK
   echo "indel realignment"
   ## get indel targets
      java \
         -Djava.io.tmpdir=${outDir}/tmp \
         -jar ${gatkDir}/GenomeAnalysisTK.jar \
         -T RealignerTargetCreator \
         -R ${ref} \
         -I ${outDir}/dedupBam/${name}.dedup.bam \
         -o ${outDir}/realignLists/${name}.list

   ## re-align around InDels
      java -Xmx20g \
         -Djava.io.tmpdir=${outDir}/tmp \
         -jar ${gatkDir}/GenomeAnalysisTK.jar \
         -T IndelRealigner \
         -R ${ref} \
         -I ${outDir}/dedupBam/${name}.dedup.bam \
         -targetIntervals ${outDir}/realignLists/${name}.list \
         -o ${finalDir}/${name}.dedup.indel.bam

