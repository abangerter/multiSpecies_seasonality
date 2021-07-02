### VISUALIZE FST RESULTS ###


### libraries
library(foreach)
library(data.table)
library(ggplot2)
library(cowplot)


### read in the data 
  ## fst info
    fst.out <- foreach(i=c("Dmel","Dsim","Dhyd","Zind","Dsuz"), .combine="rbind") %do% {
      ## read in the data -- pick which MAF 
        fst.wg <- readRDS(paste("/scratch/ab5dr/multiSpecies/analysis/poolfstat_inputDataFormat_MS_",i,"_wgFst.rds", sep=""))
      ## reformat 
        # turn into a data table 
          fst.wg.dt <- as.data.table(fst.wg[[1]], keep.rownames="pool.id.1")
          fst.wg.dt <- melt(fst.wg.dt, measure.vars=fst.wg.dt$pool.id.1, value.name="fst", variable.name="pool.id.2")
        # remove rows where ID1=ID2
          fst.wg.dt <- fst.wg.dt[pool.id.1!=pool.id.2]
        # set levels by ID1
          fst.wg.dt[,season1:= tstrsplit(pool.id.1,"_")[3]]
          fst.wg.dt[,season2:= tstrsplit(pool.id.2,"_")[3]]
          fst.wg.dt[,season1:= factor(season1, levels=c("S","F"))]
          fst.wg.dt[,season2:= factor(season2, levels=c("S","F"))]
          fst.wg.dt[,year1:= tstrsplit(pool.id.1,"_")[2]]
          fst.wg.dt[,year2:= tstrsplit(pool.id.2,"_")[2]]
          setkey(fst.wg.dt, year1, season1, year2, season2)
          fst.wg.dt[,pool.id.1:= factor(pool.id.1, levels=unique(fst.wg.dt$pool.id.1), ordered=TRUE)]
          fst.wg.dt[,pool.id.2:= factor(pool.id.2, levels=unique(fst.wg.dt$pool.id.1), ordered=TRUE)]
        # remove when ID1 < ID2 based on factor to remove doubles of comparisons (get half of the fst matrix)
          fst.wg.dt <- fst.wg.dt[! pool.id.2 < pool.id.1]
      ## add in species column and prep final table 
        fst.wg.dt <- fst.wg.dt[,c("pool.id.1", "pool.id.2", "fst")]
        fst.wg.dt[,species:=i]
      ## print and go 
        fst.wg.dt
    }
  ## time changes 
    time.dt <- fread("/scratch/ab5dr/multiSpecies/analysis/deltaT.txt", header=T)
  ## merge
    setkey(fst.out, pool.id.1, pool.id.2)
    setkey(time.dt, pool.id.1, pool.id.2)
    fst.out <- merge(fst.out, time.dt)
    
    
### prep species names for plotting 
  spec.name <- data.table(species=c("Dmel","Dsim","Dhyd","Zind","Dsuz"),
                          species.2=c("D. melanogaster","D. simulans","D. hydei","Z. indianus","D. suzukii"))
  setkey(spec.name, species)
  setkey(fst.out, species)
  fst.out <- merge(fst.out, spec.name)

  
### plot 
  ggplot(fst.out, aes(x=deltaT, y=fst, color=yearComp)) + geom_point() + theme_half_open() + facet_grid(~species.2) + 
    labs(x="Days Between Sampling", y=expression("F"["ST"]), color="Sample Year\nComparison") + panel_border()
    
  