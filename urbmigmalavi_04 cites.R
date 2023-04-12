## MalAvi, urbanization, migration
## 04_cites
## danbeck@ou.edu
## last updated 4/11/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(easyPubMed)

## load in MalAvi with traits
setwd("~/Desktop/urbmigmalavi/flat files")
data=read.csv("MalAvi hosts with traits_cleaned.csv")

## collect citations per host species, using latin AND common
cites=c()
for(i in 1:length(data$tip)) {
  
  ## split into genus and species
  x=strsplit(data$tip[i]," ")[[1]]
  
  ## paste
  x=paste(x,collapse=" AND ")
  
  ## repeat for common name
  y=strsplit(data$English[i]," ")[[1]]
  y=paste(y,collapse=" AND ")
  
  ## add parentheses
  x=paste("(",x,")",sep="")
  y=paste("(",y,")",sep="")
  
  ## combine
  xy=paste(x," AND ",y,sep="")
  
  ## get citations
  counts=as.numeric(as.character(get_pubmed_ids(xy)$Count))
  cites[i]=counts
  print(paste(data$tip[i],", cites = ",counts,", ",i,"/",nrow(data),sep=""))
}

## compile all citations
cdata=data.frame(tip=data$tip,
                 cites=cites)

## export
setwd("~/Desktop/urbmigmalavi/flat files")
write.csv(cdata,"urbmigmalavi host citations.csv")