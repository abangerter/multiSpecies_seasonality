### DATA QUALITY PLOTTING ### 


### libraries 
library(data.table)
library(gdsfmt)
library(SeqArray)
library(foreach)
library(ggplot2)
library(cowplot)


### get filtered gds files 
  for(i in c("Dmel","Dsim","Dhyd","Zind","Dsuz")) {
    print(i)
    vcf.fn <- paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,".VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt.recode.vcf",sep="")
    seqVCF2GDS(vcf.fn, paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_filtered_seqA.gds", sep=""))
  }


### Get Effective Coverage 
  ## read in and get read depth information -- loop through each species 
    rd.dt <- foreach(i=c("Dmel","Dsim","Dhyd","Zind","Dsuz"), .combine="rbind") %do% {
      print(i)
      ## open gds file 
        genofile <- seqOpen(paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_filtered_seqA.gds", sep=""))
      ## extract RD information 
        temp.rd <- seqGetData(genofile, "annotation/format/DP")
        temp.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
        colnames(temp.dt) <- seqGetData(genofile, "sample.id")
        temp.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
        temp.dt[,pos:=seqGetData(genofile, "position")]
      ## get into output long-format data table
        pops <- seqGetData(genofile, "sample.id")
        temp.dt <- melt(temp.dt, measure.vars=pops, value.name="RD", variable.name="sampleID")
      ## and go
        temp.dt
    }
  ## read in metadata 
    info.dt <- fread("/scratch/ab5dr/multiSpecies/analysis/basicInfo.txt", header=T) 
  ## merge 
    setkey(rd.dt, sampleID)
    setkey(info.dt, sampleID)
    rd.dt <- merge(rd.dt, info.dt)
  ## adjust effective coverage for known X chromsome SNPs
    rd.dt[Species=="D. melanogaster",Nchromosomes:=ifelse(chrom_scaff=="X", NfliesInPool, Nchromosomes)] 
    rd.dt[Species=="D. simulans",Nchromosomes:=ifelse(chrom_scaff=="NC_052525.1", NfliesInPool, Nchromosomes)]
    # brief aside, identifying Dsuz X scaffolds using info from S.Table 3 from Paris et al 2020
      ref.dt <- fread("/scratch/ab5dr/multiSpecies/ref/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.dict", header=F, skip=1)
      ref.dt[,length:=tstrsplit(V3,":")[2]]
      ref.dt[,scaff:=tstrsplit(V2,":")[2]]
      x_scaffold_lengths <- c("13739618","6492416","3772065","2783430","2338044","1589149","653633","90352","86498","57073","39583","33763","30305")
      ref.dt <- ref.dt[length %in% x_scaffold_lengths]
      dsuz.x <- ref.dt$scaff
    # correct Dsuz by the male vs female ratio in the sequenced samples
      rd.dt[sampleID=="MS_2017_S_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.58*NfliesInPool), Nchromosomes)] 
      rd.dt[sampleID=="MS_2017_F_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.31*NfliesInPool), Nchromosomes)] 
      rd.dt[sampleID=="MS_2018_S_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.41*NfliesInPool), Nchromosomes)] 
      rd.dt[sampleID=="MS_2018_F_Dsuz",Nchromosomes:=ifelse(chrom_scaff %in% dsuz.x, round(1.11*NfliesInPool), Nchromosomes)] 
  ## get effective coverage per SNP 
    rd.dt[,effCov:=(RD*Nchromosomes)/(RD+Nchromosomes)]
  ## get mean and stdev of effcov for each sample 
    ec.ag <- rd.dt[,list(meanEC=mean(effCov, na.rm=T),
                         stdevEC=sd(effCov, na.rm=T)),
                   by=list(sampleID, Species)]
    ec.ag[,YearCollected:=tstrsplit(sampleID,"_")[2]]


### Get summary of SNP filtering 
  snp.ag <- info.dt[,list(NSNPs_pre=mean(SNPsPreFiltering),
                          retained=mean(SNPsPostFiltering)),
                    by=list(Species)]
  snp.ag[,removed:=NSNPs_pre-retained]
  snp.ag.melt <- melt(snp.ag, measure.vars=c("retained","removed"), value.name="nSNPs", variable.name="Filtering")
  snp.ag.melt[,Filtering:= factor(Filtering, levels=c("retained","removed"))]

  
### prep for plotting 
  ec.ag[,season:=tstrsplit(sampleID,"_")[3]]
  ec.ag[,season:=ifelse(season=="S","spring","fall")]
  ec.ag[,season:=factor(season, levels=c("spring","fall"))]
  
    
### Plot Figure 2 
  ## effective coverage 
    a <- ggplot(ec.ag) + geom_point(aes(x=sampleID, y=meanEC, color=Species)) + 
      geom_linerange(aes(x=sampleID, ymin=meanEC-stdevEC, ymax=meanEC+stdevEC, color=Species)) + 
      theme_half_open() + geom_hline(yintercept=50, linetype="dashed") + ylab("Effective Coverage") + 
      facet_grid(.~YearCollected + season, scales="free_x") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
  ## snp filtering 
    b <- ggplot(snp.ag.melt, aes(x=Species, y=nSNPs, color=Filtering, fill=Filtering)) + 
      geom_bar(stat="identity", position=position_stack(reverse=TRUE)) + theme_classic() + ylab("N SNPs")
  ## plot 
    plot_grid(a,b, ncol=1, labels="AUTO")
    

    