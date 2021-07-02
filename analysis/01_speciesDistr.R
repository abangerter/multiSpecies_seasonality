### SPECIES DISTRIBUTION THROUGH TIME ###


### libraries
library(data.table)
library(ggplot2)
library(cowplot)


### read in the data 
  ## read in
    spec.dt <- fread("/scratch/ab5dr/multiSpecies/analysis/speciesCounts.txt", header=T)
  ## reformat and condense into total counts per collection
    spec.dt <- spec.dt[species!="Drosophila mel/sim"] ## unevenly counted, so remove rare counts from dataset
    spec.dt <- spec.dt[,list(TotalCount=sum(Count)),
                       by=list(Year, CollectionDate, JulianCollectionDate, species)]
  ## reformat for plotting species of interest 
    spec.vec <- c("Drosophila melanogaster", "Drosophila simulans", "Drosophila hydei","Drosophila suzukii","Zaprionus indianus")
    spec.dt[,speciesCat:=ifelse(species %in% spec.vec, species, "other")]
    spec.ag <- spec.dt[,list(TotalCount=sum(TotalCount, na.rm=T)), by=list(JulianCollectionDate, Year, speciesCat)]
    spec.ag[,date:=paste(Year, JulianCollectionDate, sep="_")]
  ## adjust melanogaster and simulans counts (only counts of males in the counts file)
    spec.ag[speciesCat %in% c("Drosophila melanogaster", "Drosophila simulans"), TotalCount:= 2*TotalCount]
  ## rename to G. species
    spec.name <- data.table(speciesCat=c("Drosophila hydei", "Drosophila melanogaster", "Drosophila simulans", "Drosophila suzukii", "Zaprionus indianus","other"),
                            species=c("D. hydei", "D. melanogaster", "D. simulans", "D. suzukii", "Z. indianus","other"))
    setkey(spec.ag, speciesCat)
    setkey(spec.name, speciesCat)
    spec.ag <- merge(spec.ag, spec.name)
  ## set levels to species
    spec.ag[,species:=factor(species, levels=c("D. hydei", "D. melanogaster", "D. simulans", "D. suzukii", "Z. indianus","other"))]


### plot with a continuous x axis if possible 
  ggplot(spec.ag, aes(x=JulianCollectionDate, y=TotalCount, color=species, fill=species)) + geom_bar(stat="identity", position="fill", width=3) + 
    facet_grid(Year~.) + theme_half_open() + labs(x="Collection Date", y="Proportion in Sample", fill="species") + guides(color=FALSE)


