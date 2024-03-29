## MalAvi, urbanization, migration
## 02_trait merge
## danbeck@ou.edu
## last updated 4/4/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(readxl)
library(reshape2)

## read in edgelist
#setwd("~/Desktop/urbmigmalavi/flat files")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
edge=read.csv("MalAvi edgelist.csv")
data=edge

## alternative lineage to NA
data$Alt_Lineage_Name=as.character(data$Alt_Lineage_Name)
data$Alt_Lineage_Name=ifelse(data$Alt_Lineage_Name=="",NA,data$Alt_Lineage_Name)

#make songbird-only edgelist?
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/AVONET")

#load in avonet data so we have orders to link up to species
avonet=read_excel("AVONET Supplementary dataset 1.xlsx",4)
taxonomy = avonet[, c("Species3", "Order3")]
taxonomy$tip = taxonomy$Species3

#merge
edge=merge(edge,taxonomy,by="tip",all.x=T)
songbird_edgelist <- edge %>% filter(Order3 == "Passeriformes")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
write.csv(songbird_edgelist,"songbirds only Malavi edgelist.csv")
rm(avonet)
rm(songbird_edgelist)

#back to merging host traits to big df 
## check name
mis=setdiff(data$tip,avonet$Species3)

## fix label
avonet$tip=avonet$Species3

## merge
data=merge(data,avonet,by="tip",all.x=T)


## collapse into host data
lset=list()
for(i in 1:length(unique(data$tip))){
  
  ## subset
  set=data[data$tip%in%unique(data$tip)[i],]
  
  ## collapse lineage
  nset=data.frame(tip=unique(data$tip)[i],
                  lineages=paste(set$Lineage_Name,collapse=", "),
                  alt_lineages=paste(set$Alt_Lineage_Name[!is.na(set$Alt_Lineage_Name)],collapse=", "))
  lset[[i]]=nset
}

## collapse
data=do.call(rbind.data.frame,lset)
rm(lset,nset,set,i)

## load AVONET
#setwd("/Users/danielbecker/OneDrive - University of Oklahoma/Becker Lab/Datasets/AVONET")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/AVONET")

avonet=read_excel("AVONET Supplementary dataset 1.xlsx",4)

## check name
mis=setdiff(data$tip,avonet$Species3)

## fix label
avonet$tip=avonet$Species3

## merge
data=merge(data,avonet,by="tip",all.x=T)

## load BirdTree taxonomy
#setwd("~/OneDrive - University of Oklahoma/Becker Lab/Datasets/BirdTree")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/BirdTree")

pnames=read.csv("BLIOCPhyloMasterTax.csv")

## check mismatch
mis=setdiff(data$tip,pnames$Scientific)

## merge
pnames$tip=pnames$Scientific
data=merge(data,pnames,by="tip",all.x=T)
rm(mis,pnames)

## Merge in Gonzalez-Lagos et al. 2022 (urban traits)
#setwd("/Users/danielbecker/OneDrive - University of Oklahoma/Becker Lab/Datasets/Gonzalez-Lagos et al 2021 Ecography")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/Gonzalez-Lagos et al 2021 Ecography")
gdata=read.csv("DataS1.csv")

## fix
gdata$tip=gsub("_"," ",gdata$animal)

## check name
mis=setdiff(data$tip,gdata$tip)

## fix synonyms
rownames(gdata)=gdata$tip

## Acrocephalus baeticatus/Acrocephalus scirpaceus
new=gdata["Acrocephalus scirpaceus",]
new$tip="Acrocephalus baeticatus"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Otus lempiji/Otus bakkamoena
new=gdata["Otus bakkamoena",]
new$tip="Otus lempiji"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Otus lettia/Otus bakkamoena
new=gdata["Otus bakkamoena",]
new$tip="Otus lettia"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Parus caeruleus/Parus teneriffae
new=gdata["Parus caeruleus",]
new$tip="Parus teneriffae"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Turdus nigriceps/Turdus subalaris
new=gdata["Turdus nigriceps",]
new$tip="Turdus subalaris"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Poospiza nigrorufa/Poospiza whitii
new=gdata["Poospiza nigrorufa",]
new$tip="Poospiza whitii"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Momotus aequatorialis/Momotus momota
new=gdata["Momotus momota",]
new$tip="Momotus aequatorialis"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Nectarinia ludovicensis/Cinnyris ludovicensis
new=gdata["Cinnyris ludovicensis",]
new$tip="Nectarinia ludovicensis"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## Phylloscopus nitidus/Phylloscopus trochiloides
new=gdata["Phylloscopus trochiloides",]
new$tip="Phylloscopus nitidus"
rownames(new)=new$tip
gdata=rbind.data.frame(gdata,new)

## check name
mis=setdiff(data$tip,gdata$tip)
rm(mis,new)

## clean trait data
gdata=gdata[c("tip","urban","humanDisturbed")]

## merge
data=merge(data,gdata,by="tip",all.x=T)
rm(gdata,avonet)

## as factor
data$urban=factor(data$urban)
data$humanDisturbed=factor(data$humanDisturbed)
data$Migration=factor(data$Migration)

##trim to Passeriformes
library(dplyr)
passer_data <- data %>% filter(Order3 == "Passeriformes")


## export cleaned data
#setwd("~/Desktop/urbmigmalavi/flat files")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
write.csv(data,"MalAvi hosts with traits_cleaned_Passeriformes.csv")

## pull lineages
ldata=edge
ldata=ldata[!duplicated(ldata$Lineage_Name),]
ldata=ldata[c("Lineage_Name","Alt_Lineage_Name","parasiteGenus")]

#same for trimmed edgelist
tldata=songbird_edgelist
tldata=tldata[!duplicated(tldata$Lineage_Name),]
tldata=tldata[c("Lineage_Name","Alt_Lineage_Name","parasiteGenus")]

## export
setwd("~/Desktop/urbmigmalavi/flat files")
write.csv(ldata,"MalAvi lineages cleaned.csv")
write.csv(tldata,"songbirds only Malavi lineages cleaned.csv")
