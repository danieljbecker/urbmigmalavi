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
library(paco)
library(Rtapas)

## load edgelist
setwd("~/Desktop/urbmigmalavi/flat files")
edge=read.csv("MalAvi edgelist.csv")

## for trialing out code, subsample the egdelist
set.seed(1)
edge=edge[sample(nrow(edge),0.25*nrow(edge)),]

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
ptree=keep.tip(ptree,spec)
rm(spec)

## make a H-P association matrix
hptab=with(edge,table(tiplab,Lineage_Name))
hptab=table(edge$tiplab,edge$Lineage_Name)

## binary
hptab=ifelse(hptab>0,1,0)

# melt and visualize the matrix (computationally taxing)
hptab_melt=melt(hptab)
names(hptab_melt)=c("tiplab","Lineage_Name","value")
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

## check
length(colnames(htree_dist))==length(unique(edge$tiplab))
length(colnames(ptree_dist))==length(unique(edge$Lineage_Name))

# ## Rtapas setup
# N=100
# nset=one2one_f(hptab,reps=N,interval=c(1,10),plot=T)
# pmax=max_cong(HS=hptab,
#               treeH=htree,
#               treeS=ptree,
#               n=nset,
#               N=N,
#               method="paco")
# 
# ## simple tangelgram
# tangle_gram(treeH=htree,
#             treeS=ptree,
#             HS=hptab,
#             fqtab=pmax,
#             colscale="diverging",
#             colgrad=viridis(10))

## parafit
# set.seed(1)
# pfit=parafit(host.D=htree_dist,para.D=ptree_dist,HP=hptab,
#              correction="cailliez",nperm=10,test.links=TRUE)

## paco
D=prepare_paco_data(H=htree_dist,P=ptree_dist,HP=hptab)
D=add_pcoord(D,correction="cailliez")
pac=PACo(D,nperm=999,seed=1,method="r0",symmetric=F) ## assumes column group tracks row group

## get interaction-specific cophylogenetic contributions based on jacknife
pac_links=paco_links(pac)

## pull out contributions
res=residuals_paco(pac_links$proc)
rdata=data.frame(res)
rdata$pairs=rownames(rdata)
x=strsplit(rdata$pairs,"-")
rdata$tiplab=sapply(x,function(x) x[1])
rdata$Lineage_Name=sapply(x,function(x) x[2])

## jackknife
jdata=data.frame(pac_links$jackknife)
names(jdata)=c("jackknife")
jdata$pairs=rownames(jdata)

## combine
rdata=merge(rdata,jdata,by="pairs",all.x=T)
rm(jdata)

## weight as rescaled res
wei=plotrix::rescale(res,c(2,0.1)) ## large res = little weight
wei=((res^-2)/4)

## make the interaction matrix for cophyloplot format
imat=hptab_melt
imat=imat[which(imat$value>0),]
imat$value=NULL

## cophyloplot
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
cophyloplot(htree,ptree,assoc=imat,show.tip.label=F,
            use.edge.length=F,lwd=0.01,space=2000,gap=5,length.line=-5)
## lwd = wei