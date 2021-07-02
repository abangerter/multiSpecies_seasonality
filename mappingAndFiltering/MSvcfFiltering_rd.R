#! /bin/usr/env Rscript
### FILTERING VCF FILES FOR READ DEPTH OUTLIERS ###


### libraries
library(gdsfmt)
library(SeqArray)
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)


### make gds files
  for(i in c("Dmel","Dsim","Dhyd","Zind","Dsuz")) {
    print(i)
    vcf.fn <- paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,".VarScan.newHeader.noIndelReg.noRep.biAllelic.recode.vcf",sep="")
    seqVCF2GDS(vcf.fn, paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_preRDfilt_seqA.gds", sep=""))
  }


### choosing filter threshold 
  ## loop
    rd.out <- foreach(i=c("Dmel","Dsim","Dhyd","Zind","Dsuz"), .combine="rbind") %do% {
      print(i)
      ## read in gds file 
        genofile <- seqOpen(paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_preRDfilt_seqA.gds", sep=""))
      ## pull out read depths 
        # total RD
          temp.rd <- seqGetData(genofile, "annotation/format/DP")
          RD.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
          colnames(RD.dt) <- seqGetData(genofile, "sample.id")
          RD.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
          RD.dt[,pos:=seqGetData(genofile, "position")]
         # melt table
          pops <- seqGetData(genofile, "sample.id")
          RD.dt <- melt(RD.dt, measure.vars=pops, value.name="RD", variable.name="pop.id")
        # aggregate site reads depths 
          RD.ag <- RD.dt[,list(totalRD=sum(RD, na.rm=T)),
                               by=list(chrom_scaff, pos)]
      ## close genofile
        seqClose(genofile)
      ## out and next
        RD.ag[,species:=i]
        RD.ag
    }
  ## plot RD histograms 
    ggplot(rd.out, aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free")
  ## remove Dmel and Dsim not autosomes (2L,2R,3L,3R) and X
    rd.sub <- rbindlist(list(rd.out[species=="Dmel"][chrom_scaff %in% c("X","2L","2R","3L","3R")],
                             rd.out[species=="Dsim"][chrom_scaff %in% c("NC_052520.1","NC_052521.1","NC_052522.1","NC_052523.1","NC_052525.1")],
                             rd.out[species=="Dhyd"], rd.out[species=="Zind"], rd.out[species=="Dsuz"]))
  ## re-plot RD histograms to see if removing those non chromosomes scaffolds dramatically shift RD distr or not
    ggplot(rd.sub, aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free")
  ## pick a few thresholds based on total RD
    thresholds <- rd.sub[,list(threshold.99=quantile(totalRD, 0.99, na.rm=T),
                               threshold.95=quantile(totalRD, 0.95, na.rm=T),
                               threshold.90=quantile(totalRD, 0.90, na.rm=T)),
                         by=list(species)]
  ## plot -- "zoomed in" to see the main distribution and where the thresholds land relative to the bulk of the data 
    a <- ggplot(rd.sub[species=="Dhyd"][totalRD<1000], aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free") +
      geom_vline(xintercept=thresholds[species=="Dhyd"]$threshold.99, color="red") + 
      geom_vline(xintercept=thresholds[species=="Dhyd"]$threshold.95, color="blue") + 
      geom_vline(xintercept=thresholds[species=="Dhyd"]$threshold.90, color="green")
    b <- ggplot(rd.sub[species=="Dmel"][totalRD<1000], aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free") +
      geom_vline(xintercept=thresholds[species=="Dmel"]$threshold.99, color="red") + 
      geom_vline(xintercept=thresholds[species=="Dmel"]$threshold.95, color="blue") + 
      geom_vline(xintercept=thresholds[species=="Dmel"]$threshold.90, color="green")
    c <- ggplot(rd.sub[species=="Dsim"][totalRD<1000], aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free") +
      geom_vline(xintercept=thresholds[species=="Dsim"]$threshold.99, color="red") + 
      geom_vline(xintercept=thresholds[species=="Dsim"]$threshold.95, color="blue") + 
      geom_vline(xintercept=thresholds[species=="Dsim"]$threshold.90, color="green")
    d <- ggplot(rd.sub[species=="Dsuz"][totalRD<1200], aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free") +
      geom_vline(xintercept=thresholds[species=="Dsuz"]$threshold.99, color="red") + 
      geom_vline(xintercept=thresholds[species=="Dsuz"]$threshold.95, color="blue") + 
      geom_vline(xintercept=thresholds[species=="Dsuz"]$threshold.90, color="green")
    e <- ggplot(rd.sub[species=="Zind"][totalRD<1000], aes(x=totalRD)) + geom_histogram() + facet_wrap(~species, scales="free") +
      geom_vline(xintercept=thresholds[species=="Zind"]$threshold.99, color="red") + 
      geom_vline(xintercept=thresholds[species=="Zind"]$threshold.95, color="blue") + 
      geom_vline(xintercept=thresholds[species=="Zind"]$threshold.90, color="green")
    plot_grid(a,b,c,d,e, labels=NULL)
  ### selecting 99% threshold ###


### RD based filtering 
  for(i in c("Dmel","Dsim","Dhyd","Zind","Dsuz")) {
    print(i)
    ## read in gds file 
      genofile <- seqOpen(paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_preRDfilt_seqA.gds", sep=""))
    ## pull out read depths 
      # total RD
        temp.rd <- seqGetData(genofile, "annotation/format/DP")
        RD.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
        colnames(RD.dt) <- seqGetData(genofile, "sample.id")
        RD.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
        RD.dt[,pos:=seqGetData(genofile, "position")]
      # melt table
        pops <- seqGetData(genofile, "sample.id")
        RD.dt <- melt(RD.dt, measure.vars=pops, value.name="RD", variable.name="pop.id")
    ## do some if/else trickery to match filtering above to determine thresholds 
      if(i=="Dmel") {
        RD.dt <- RD.dt[chrom_scaff %in% c("X","2L","2R","3L","3R")]
      } else if(i=="Dsim") {
        RD.dt <- RD.dt[chrom_scaff %in% c("NC_052520.1","NC_052521.1","NC_052522.1","NC_052523.1","NC_052525.1")]
      } else {
        print("no additional filtering needed")
      }
    ## identify the quantile thresholds for filtering 
      totalRD <- RD.dt[,list(totalRD=sum(RD, na.rm=T)),
                       by=list(chrom_scaff, pos)]
      threshold <- quantile(totalRD$totalRD, 0.99)
      totalRD <- totalRD[totalRD <= threshold][,c("chrom_scaff","pos")]
    ## filter table be quantile threshold 
      setkey(totalRD, chrom_scaff, pos)
      setkey(RD.dt, chrom_scaff, pos)
      RD.dt <- merge(RD.dt, totalRD)
    ## write out chromosome/scaffold-position information for SNPs to keep in filtering
      write.table(RD.dt[,c("chrom_scaff", "pos")], paste("/scratch/ab5dr/multiSpecies/vcfs/MS_",i,"_RDkeep.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
    ## close out gds file 
      seqClose(genofile)
  }




    
    
    
######## deleted code #####
    # pull out alternative allele read depth 
    ad <- seqGetData(genofile, var.name="annotation/format/AD", .useraw=FALSE)
    ad.dt <- foreach(j=1:dim(ad$data)[1], .combine="cbind") %do% {data.table(ad$data[j,])}
    colnames(ad.dt) <- seqGetData(genofile, "sample.id")
    ad.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
    ad.dt[,pos:=seqGetData(genofile, "position")]
    # pull out reference allele read depth 
    ref <- seqGetData(genofile, var.name="annotation/format/RD", .useraw=FALSE)
    ref.dt <- foreach(j=1:dim(ref$data)[1], .combine="cbind") %do% {data.table(ref$data[j,])}
    colnames(ref.dt) <- seqGetData(genofile, "sample.id")
    ref.dt[,chrom_scaff:=seqGetData(genofile, "chromosome")]
    ref.dt[,pos:=seqGetData(genofile, "position")]
    # melt tables 
    pops <- seqGetData(genofile, "sample.id")
    RD.dt <- melt(RD.dt, measure.vars=pops, value.name="RD", variable.name="pop.id")
    ad.dt <- melt(ad.dt, measure.vars=pops, value.name="NaltReads", variable.name="pop.id")
    ref.dt <- melt(ref.dt, measure.vars=pops, value.name="NrefReads", variable.name="pop.id")
    
