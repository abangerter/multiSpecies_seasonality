### DRAFT README 


1. SPECIES DISTRIBUTION - FIGURE 1
	The script 01_speciesDistr.R uses metadata file speciesCounts.txt to generate figure 1.

2. DATA QUALITY - FIGURE 2
	The script 02_dataQsumm.R generates GDS files from the filtered VCF file, calculates effective coverage, summarizes filtering, and generates figure 2. It uses metadata file basicInfo.txt 

3. FST - FIGURE 3
	A. The script 03A_fst.R estimates whole-genome, pairwise FST within each species using the filtering VCF file. 
	B. The script 03B_fst_plotting.R takes the FST estimates generated in step A and generates Figure 3. It uses metadata file deltaT.txt

4. LINEAR MODEL - FIGURE 4
	A. The script 04A_lm_dataPrep.R formats and prepares data for running through the linear model. This script uses gds files generated by 02_dataQsumm.R and metadata file basicInfo.txt
	B. The script 04B_lm_forRivanna.R runs the linear model of season on allele frequency changes per SNP per species. In it, it contains the options to run on both the regular data and a highly filtered dataset. 
	C. The script 04C_lm_outputAnalysis.R takes the output of the linear model, does p-value correction, generates Table 3 and Figure 4. It also contains options for both the regular data and the highly filtered dataset, though the highly filtered dataset is what was used in the manuscript. 

