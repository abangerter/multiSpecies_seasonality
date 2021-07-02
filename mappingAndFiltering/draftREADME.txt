#### DRAFT README


1. MAPPING RAW READS
	A. First run map_reads.sh. This is accompanied by MSlist.txt which lists all the sample IDs and map_reads.slurm, which combined is how we submitted these scripts on our computing cluster. 
	B. Second, run map_reads_merge.sh. This scripts has annotations within about how to modify it based on if it's being run on samples that were sequenced in two different sequencing runs. Only a subset of samples, found in MSlist_merge.txt, needed to be merged as only a subset of samples were sequenced on two different sequencing runs. Those samples only sequenced in one sequencing run can be found in MSlist_noMerge.txt. To submit these scripts, we used map_reads_merge.slurm with the two different input lists of samples in there that just need to be un/commented out to swap out what list is being submitted. 

2. DECONTAMINATE MAPPED READS
	Run ms_decontaminate.sh - but there is a large chunk of code that should only be run once to build up the combined "reference" genomes to be used in the competitive mapping. Also, you will need the python scripts in the pythonScripts directory. This is run using the MSlist_decontaminate.txt and ms_process.slurm to work on our computing cluster. 

3. CALL VARIANTS
	Variants are called using MSvarScan.sh, which includes both a samtools mpileup step and the varscan variant calling. You will also need the "bam lists" to make sure you know the input order into the varscan for the reheader step in 4B. It should be submitted with the species "code" as an argument:
		./MSvarScan.sh Dmel
		./MSvarScan.sh Dsim
		./MSvarScan.sh Dsuz
		./MSvarScan.sh Dhyd
		./MSvarScan.sh Zind

4. FILTER VARIANTS 
	A) First, run MSvarScan_indels.sh to obtain the deNovo indel calls. It should be submitted similarly to the above step:
		./MSvarScan_indels.sh Dmel
		./MSvarScan_indels.sh Dsim
		./MSvarScan_indels.sh Dsuz
		./MSvarScan_indels.sh Dhyd
		./MSvarScan_indels.sh Zind
	B) Second, run reheader.sh to rename the sample header. VarScan automatically names samples in the order from the bamlist input file. The reheader files in the reheaderTextFiles directory make sure that SampleN are matched properly. 
	C) Third, run MSvcfFiltering.sh. Again, run it like above scripts with the species code as an argument. Keep in mind that the bed file generated in the first step can result in negative positions. This was fixed by hand in a text editor once the error was identified. The repeat masker bed files are found in the repeatMaskerFiles directory. 
	D) Fourth, run MSvcfFiltering_rd.R to obtain the list of SNPs to keep based on read depth filtering. 
	D) Fifth, run msRDfiltering.sh to get the final filtered VCF file. 



