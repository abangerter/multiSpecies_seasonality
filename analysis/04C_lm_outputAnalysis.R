### ANALYSIS OF LM OUTPUT TO IDENTIFY SEASONAL SNPS ###


### libraries 
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)


### read in the data & get p-value corrections 
  ## MAF 0.01
    # lm.out <- foreach(i=c("Dmel","Dsim","Dhyd","Zind","Dsuz"), .combine="rbind") %do% {
    #   temp <- readRDS(paste("/scratch/ab5dr/multiSpecies/analysis/lm_out_",i,".rds", sep=""))
    #   temp[,fdr:=p.adjust(p, method="fdr")]
    #   temp[,bonf_p:=p.adjust(p, method="bonferroni")]
    #   temp
    # }
  ## MAF 0.1
    lm.out <- foreach(i=c("Dmel","Dsim","Dhyd","Zind","Dsuz"), .combine="rbind") %do% {
      temp <- readRDS(paste("/scratch/ab5dr/multiSpecies/analysis/lm_out_highFilter_",i,".rds", sep=""))
      temp[,fdr:=p.adjust(p, method="fdr")]
      temp[,bonf_p:=p.adjust(p, method="bonferroni")]
      temp
    }


### look at distribution of p values pre and post adjustment
  ggplot(lm.out, aes(x=p)) + geom_histogram() + facet_wrap(~species)
  ggplot(lm.out, aes(x=fdr)) + geom_histogram() + facet_wrap(~species)
  ggplot(lm.out, aes(x=bonf_p)) + geom_histogram() + facet_wrap(~species)
  
  
### how many snps with p <0.05?
  lm.ag <- lm.out[,list(Ntotal=.N,
                        Np=sum(p<0.05, na.rm=T),
                        Nfdr0.05=sum(fdr<0.05, na.rm=T),
                        Nfdr0.1=sum(fdr<0.1, na.rm=T),
                        Nbonf=sum(bonf_p<0.05, na.rm=T)),
                  by=list(species)]
  ## MAF 0.01
    #            species     Np Nfdr0.05 Nfdr0.1 Nbonf
    # 1: D. melanogaster  81666        0       0     0
    # 2:     D. simulans 269519      218     960     8
    # 3:        D. hydei 101070        5       9     5
    # 4:     Z. indianus 323626       28      52     6
    # 5:      D. suzukii 302161        5       5     4
  ## MAF 0.1 -- TABLE 3
    #            species  Ntotal     Np Nfdr0.05 Nfdr0.1 Nbonf
    # 1: D. melanogaster 1169311  72040        0       0     0
    # 2:     D. simulans 2345168 247340      415    2837     6
    # 3:        D. hydei 1510790  88323        8      10     5
    # 4:     Z. indianus 3588541 287021       29     111     6
    # 5:      D. suzukii 4528392 224524        2       2     1
  
  
### try to look at AF changes of "significant" SNPs
  ## read back in AF data 
    af.dt <- readRDS("/scratch/ab5dr/multiSpecies/analysis/rawAFdata_forLM.rds")
    af.dt[,snpID:=paste(chrom_scaff, pos, sep="_")]
    af.dt <- af.dt[,c("sampleID","snpID","spring_fall","YearCollected","af_cor")]
  ## top 500 SNPs per species 
    lm.sig <- foreach(i=c("D. melanogaster","D. simulans","D. hydei","Z. indianus","D. suzukii"), .combine="rbind") %do% {
      # subset 
        temp <- lm.out[species==i]
      # merge with AF data 
        setkey(temp, snpID)
        setkey(af.dt, snpID)
        temp <- merge(temp, af.dt)
      # get rid of SNPs where AF hits 1 or 0
        snpIDs <- c(unique(temp[af_cor==0]$snpID), unique(temp[af_cor==1]$snpID))
        temp <- temp[! snpID %in% snpIDs]
      # grab the top 500 SNPs
        setkey(temp, p)
        n <- length(unique(temp$sampleID))
        temp.out <- head(temp, 500*n)
      # out and onto the next
        temp.out
    }
  ## merge information 
    lm.sig[,sampleTimePoint:=paste(spring_fall, YearCollected, sep="_")]
    lm.sig[,sampleTimePoint:=factor(sampleTimePoint, levels=c("spring_2017","fall_2017","spring_2018","fall_2018"))]
  ## use beta to polarize allele frequency
    lm.sig[,af.adj:=ifelse(beta>0, af_cor, 1-af_cor)]
  ## plot S/F allele frequencies of "significant" SNPs -- FIGURE 4
     ggplot(lm.sig, aes(x=sampleTimePoint, y=af.adj, group=snpID)) + geom_line(alpha=0.2) + facet_wrap(~species) + theme_half_open() +
      scale_x_discrete(labels=c("Spring \n2017","Fall \n2017","Spring \n2018","Fall \n2018")) + panel_border() +
      labs(x="Collection Period", y="Allele Frequency")
    

  