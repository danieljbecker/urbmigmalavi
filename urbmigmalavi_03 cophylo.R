## MalAvi, urbanization, migration
## 03_cophylo
## danbeck@ou.edu
## last updated 3/26/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(ape)
library(caper)
library(reshape2)
library(ggplot2)

## load edgelist
setwd("~/Desktop/urbmigmalavi/flat files")
edge=read.csv("MalAvi edgelist.csv")

## load BirdTree
setwd("/Users/danielbecker/Desktop/urbmigmalavi/BirdTree")
htree=readRDS("BirdTree 2K tree consensus.rds")

## make tip label
edge$tiplab=gsub(" ","_",edge$tip)

## trim to species in edgelist
setdiff(edge$tiplab,htree$tip.label)
htree=keep.tip(htree,unique(edge$tiplab))

## load MalAvi phylogeny
setwd("/Users/danielbecker/Desktop/urbmigmalavi/MalAvi phylo")
ptree=read.tree("FastTree_output_tree.nhx")

## check correct match
spec=unique(edge$Lineage_Name)
setdiff(spec,ptree$tip.label)
rm(spec)

## make a H-P association matrix
hptab=with(edge,table(tiplab,Lineage_Name))

## binary
hptab=ifelse(hptab>0,1,0)

## melt and visualize the matrix (computationally taxing)
# hptab_melt=melt(hptab)
# ggplot(hptab_melt,aes(tiplab,Lineage_Name))+
#   geom_tile(aes(fill=value))+
#   scale_fill_gradient(low="wheat1",high="steelblue")+
#   theme_bw()+
#   coord_flip()+
#   theme(axis.text=element_blank(),
#         axis.title=element_blank())+
#   guides(fill=F)

## phylogenetic distances for host and parasite
htree_dist=cophenetic.phylo(htree)
ptree_dist=cophenetic.phylo(ptree)

## parafit
pfit=parafit(host.D=htree_dist,para.D=ptree_dist,HP=hptab,
             correction="cailliez",nperm=999,test.links=TRUE)
pfit