### LINEAR MODEL FOR SEASONALITY ###


### libraries 
library(gdsfmt)
library(SeqArray)
library(data.table)
library(foreach)



### pull out read depth for samples and information about number of chromosomes and build a table 
  ## loop through to get read depth information
    rd.dt <- foreach(i=c("Dmel","Dsim","Dhyd","Zind","Dsuz"), .combine="rbind") %do% {
      print(i)
      # open gds file 
        genofile <- seqOpen(paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_filtered_seqA.gds", sep=""))
      # extract total RD information 
        temp.rd <- seqGetData(genofile, "annotation/format/DP")
        temp.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
        colnames(temp.dt) <- seqGetData(genofile, "sample.id")
        temp.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
        temp.dt[,pos:=seqGetData(genofile, "position")]
      # extract alternate allele RD information 
        ad <- seqGetData(genofile, var.name="annotation/format/AD", .useraw=FALSE)
        ad.dt <- foreach(j=1:dim(ad$data)[1], .combine="cbind") %do% {data.table(ad$data[j,])}
        colnames(ad.dt) <- seqGetData(genofile, "sample.id")
        ad.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
        ad.dt[,pos:=seqGetData(genofile, "position")]
      # get into output long-format data table
        pops <- seqGetData(genofile, "sample.id")
        temp.dt <- melt(temp.dt, measure.vars=pops, value.name="RD", variable.name="sampleID")
        ad.dt <- melt(ad.dt, measure.vars=pops, value.name="NaltReads", variable.name="sampleID")
        setkey(temp.dt, sampleID, chrom_scaff, pos)
        setkey(ad.dt, sampleID, chrom_scaff, pos)
        temp.dt <- merge(temp.dt, ad.dt)
      # and go
        temp.dt
    }
  ## read in metadata 
    info.dt <- fread("/scratch/ab5dr/multiSpecies/analysis/basicInfo.txt", header=T) 
  ## merge 
    setkey(rd.dt, sampleID)
    setkey(info.dt, sampleID)
    af.dt <- merge(rd.dt, info.dt)
  ## subset down to important columns only 
    af.dt <- af.dt[,c("sampleID","chrom_scaff","pos","RD","NaltReads","Species","spring_fall","YearCollected","NfliesInPool","Nchromosomes")]
  ## write out this table
    saveRDS(af.dt, "/scratch/ab5dr/multiSpecies/analysis/RDdata.rds")


### get effective coverage corrected allele frequencies 
  af.dt <- readRDS("/scratch/ab5dr/multiSpecies/analysis/RDdata.rds")
  af.dt[,af:=NaltReads/RD]
  af.dt[Species=="Dmel",Nchromosomes:=ifelse(chrom_scaff=="X", NfliesInPool, Nchromosomes)] 
  af.dt[Species=="Dsim",Nchromosomes:=ifelse(chrom_scaff=="NC_052525.1", NfliesInPool, Nchromosomes)]
  ## brief aside, identifying Dsuz X scaffolds using info from S.Table 3 from Paris et al 2020
    ref.dt <- fread("/scratch/ab5dr/multiSpecies/ref/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.dict", header=F, skip=1)
    ref.dt[,length:=tstrsplit(V3,":")[2]]
    ref.dt[,scaff:=tstrsplit(V2,":")[2]]
    x_scaffold_lengths <- c("13739618","6492416","3772065","2783430","2338044","1589149","653633","90352","86498","57073","39583","33763","30305")
    ref.dt <- ref.dt[length %in% x_scaffold_lengths]
    dsuz.x <- ref.dt$scaff
  af.dt[sampleID=="MS_2017_S_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.58*NfliesInPool), Nchromosomes)] 
  af.dt[sampleID=="MS_2017_F_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.31*NfliesInPool), Nchromosomes)] 
  af.dt[sampleID=="MS_2018_S_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.41*NfliesInPool), Nchromosomes)] 
  af.dt[sampleID=="MS_2018_F_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.11*NfliesInPool), Nchromosomes)] 
  af.dt[,dp_cor:= floor((RD*Nchromosomes)/(RD+Nchromosomes))]
  af.dt[dp_cor<0,dp_cor:= 0]
  af.dt[,af_cor:= round(af*dp_cor)/dp_cor]
  af.dt[is.na(af_cor),af_cor:= 0]

  
### do a quick filter of corrected allele frequencies & write out files for input data into model 
  af.dt[,any_af_fixed:= any(af_cor<0.01 | af_cor>0.99), by=list(chrom_scaff, pos)] 
  af.dt[,mean_af_extreme:= (mean(af_cor)<0.10 | mean(af_cor)>0.90), by=list(chrom_scaff, pos)]
  saveRDS(af.dt, "/scratch/ab5dr/multiSpecies/analysis/rawAFdata_forLM.rds")
  af.filt.dt <- af.dt[!any_af_fixed & !mean_af_extreme]
  saveRDS(af.filt.dt, "/scratch/ab5dr/multiSpecies/analysis/filteredAFdata_forLM.rds")
  



