### RUN FST ON EACH SPECIES ### 


### library
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(poolfstat)
library(data.table)


### prep
  ## pool size vectors
    poolsize.Dmel <- c(598, 162, 240)
    poolsize.Dsim <- c(606, 130, 164)
    poolsize.Dhyd <- c(42, 128, 100)
    poolsize.Dsuz <- c(92, 72, 74, 64)
    poolsize.Zind <- c(744, 378, 252)
  ## pool name vectors
    poolnames.Dmel <- c("MS_2017_F_Dmel","MS_2018_F_Dmel","MS_2017_S_Dmel")
    poolnames.Dsim <- c("MS_2018_F_Dsim","MS_2017_F_Dsim","MS_2018_S_Dsim")
    poolnames.Dhyd <- c("MS_2017_F_Dhyd","MS_2018_S_Dhyd","MS_2018_F_Dhyd")
    poolnames.Dsuz <- c("MS_2018_F_Dsuz","MS_2017_S_Dsuz","MS_2017_F_Dsuz","MS_2018_S_Dsuz")
    poolnames.Zind <- c("MS_2018_F_Zind","MS_2017_S_Zind","MS_2017_F_Zind")
  ## put into table 
    pool.dt <- data.table(poolSize=c(poolsize.Dmel, poolsize.Dsim, poolsize.Dhyd, poolsize.Dsuz, poolsize.Zind),
                          poolNames=c(poolnames.Dmel, poolnames.Dsim, poolnames.Dhyd, poolnames.Dsuz, poolnames.Zind),
                          species=rep(c("Dmel", "Dsim", "Dhyd", "Dsuz", "Zind"), times=c(3,3,3,4,3)))

  
### loop through calculations
  for(i in c("Dmel","Dsim","Dhyd","Zind","Dsuz")) {
    print(i)
    ## select proper vectors
      poolsize <- pool.dt[species==i]$poolSize
      poolname <- pool.dt[species==i]$poolNames
    ## get data formatted
      pooldata <- vcf2pooldata(vcf.file=paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,".VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt.recode.vcf",sep=""),
                               poolsizes=poolsize, poolnames=poolname, min.maf=0.005)
      saveRDS(pooldata, file=paste("/scratch/ab5dr/multiSpecies/analysis/poolfstat_inputDataFormat_MS_",i,".rds",sep=""))
    ## calculate pairwise whole genome FST 
      pooldata <- readRDS(paste("/scratch/ab5dr/multiSpecies/analysis/poolfstat_inputDataFormat_MS_",i,".rds", sep=""))
      fst.wg <- computePairwiseFSTmatrix(pooldata, min.maf=0.005)
      saveRDS(fst.wg, file=paste("/scratch/ab5dr/multiSpecies/analysis/poolfstat_inputDataFormat_MS_",i,"_wgFst.rds", sep=""))
  }
   
   
      
      