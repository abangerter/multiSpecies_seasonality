#! /usr/bin/env Rscript

### CONTAINS TWO OPTIONS FOR USING RAW DATA VS FILTERED DATA 

### libraries
library(data.table)
library(doMC)
library(foreach)


### bring in species argument
  args <- commandArgs(trailingOnly=TRUE)
  species <- args[1]
  spec.dt <- data.table(argSpec=c("Dmel", "Dsim", "Dhyd", "Zind", "Dsuz"),
                        nameForLoop=c("D. melanogaster","D. simulans","D. hydei","Z. indianus","D. suzukii"))
  i <- spec.dt[argSpec==species]$nameForLoop

  
### prep parallel 
  registerDoMC(10)

  
### read back in data 
  #af.dt <- readRDS("/scratch/ab5dr/multiSpecies/analysis/rawAFdata_forLM.rds")
  af.dt <- readRDS("/scratch/ab5dr/multiSpecies/analysis/filteredAFdata_forLM.rds")
  
### subset down to species for iteration
  af.dt <- af.dt[Species==i]
  
### assign seasonal variable 
  af.dt[,seasonVar:= ifelse(spring_fall=="spring", 0, 1)]
  
### get rid of SNPs with missing data 
  af.dt <- af.dt[! is.na(RD)]
  
### get snpID
  af.dt[,snpID:= paste(chrom_scaff,pos, sep="_")]
  
### remove when nsamps < 3
  af.ag <- af.dt[,list(.N), by=(snpID)]
  tooFew <- af.ag[N<3]$snpID
  af.dt <- af.dt[! snpID %in% tooFew]
  
### run the linear model
  snpIDs <- unique(af.dt$snpID)
  print("starting LM")
  rd.dt <- foreach(j=snpIDs, .combine="rbind") %dopar% {
    nsamps <- length(af.dt[snpID==j]$sampleID)
    ## run the glm
      temp.glm <- summary(glm(af_cor~seasonVar, data=af.dt[snpID==j], family="binomial", weights=dp_cor))
    ## get output table ready
      temp.out <- data.table(species=i,
                             snpID=j,
                             beta=temp.glm$coefficients[2,1],
                             se=temp.glm$coefficients[2,2],
                             t=temp.glm$coefficients[2,3],
                             p=temp.glm$coefficients[2,4])
    ## go and next
      temp.out
  }
### write out the output
  print("writing out table")
  #saveRDS(rd.dt, paste("/scratch/ab5dr/multiSpecies/analysis/lm_out_",species,".rds", sep=""))
  saveRDS(rd.dt, paste("/scratch/ab5dr/multiSpecies/analysis/lm_out_highFilter_",species,".rds", sep=""))
  
