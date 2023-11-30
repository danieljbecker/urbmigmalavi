## MalAvi, urbanization, migration
## 05_pcms
## danbeck@ou.edu
## last updated 4/11/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(ape)
library(caper)
library(reshape2)
library(ggplot2)
library(brms)
library(bayestestR)
library(lme4)

## load in MalAvi with traits
#setwd("~/Desktop/urbmigmalavi/flat files")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
data=read.csv("MalAvi hosts with traits_cleaned_Passerellidae.csv")

## load in citations
cdata=read.csv("urbmigmalavi host citations.csv")
cdata$X=NULL

## merge
data=merge(data,cdata,by="tip")
rm(cdata)

## log cites
data$lcites=log1p(data$cites)

## get lineage richness
data$richness=sapply(strsplit(data$lineages,", "),function(x) length(unique(x)))

## urban, human, and migration as factors
data$urban=factor(data$urban)
data$humanDisturbed=factor(data$humanDisturbed)
data$Migration=factor(data$Migration)

## load BirdTree
#setwd("/Users/danielbecker/Desktop/urbmigmalavi/BirdTree")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/BirdTree")
htree=readRDS("BirdTree 2K tree consensus.rds")

## make tip label
data$tiplab=data$TipLabel

## trim to species in edgelist
setdiff(data$tiplab,htree$tip.label)
htree=keep.tip(htree,unique(data$tiplab))

## make label
data$label=data$tiplab

## observation
data$obs=factor(1:nrow(data))

## merge
cdata=comparative.data(phy=htree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## relabel
cdata$data$label=cdata$data$tiplab

## vis richness
ggplot(cdata$data[!is.na(cdata$data$urban),],
       aes(urban,richness,fill=Migration))+geom_boxplot()+scale_y_log10()

#richness in human-disturbed
ggplot(cdata$data[!is.na(cdata$data$humanDisturbed),],
       aes(humanDisturbed,richness,fill=Migration))+geom_boxplot()+scale_y_log10()

## compute correlation matrix
cmatrix=vcv.phylo(cdata$phy,cor=T)

## disable
options(buildtools.check = function(action) TRUE)

## fit model
set.seed(1)
bmod=brm(richness~humanDisturbed+lcites+
           (1|obs)+(1|gr(label,cov=A)),
         chain=1,iter=3000,thin=1,warmup=1500,
         data=cdata$data,family=poisson,
         data2=list(A=cmatrix))

## save HDI
hs=data.frame(bayestestR::hdi(bmod,ci=0.95))
summary(bmod)
plot(conditional_effects(bmod), points = TRUE)

## glmer (poisson with observation-level RE)
#mod=glmer(richness~urban*Migration+lcites+(1|obs),data=data,family=poisson)
#Anova(mod)

## glm.nb
#mod=glm.nb(richness~urban*Migration+lcites,data=data)
#Anova(mod)

## load in paco results
#setwd("~/Desktop/urbmigmalavi/flat files")
setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
rdata=read.csv("PACo results_Passerellidae.csv")
rdata$X=NULL

## merge with traits
rdata=merge(rdata,data,by="tiplab",all.x=T)

## vis jacknife
ggplot(rdata[!is.na(rdata$humanDisturbed),],
       aes(humanDisturbed,jackknife))+geom_boxplot()

hist(log(rdata$jackknife))

## lmer: residuals are individual links, not averaged by host sp.
mod=lmer(log(jackknife)~humanDisturbed+(1|Lineage_Name),
          data=rdata)
summary(mod)
#Anova(mod)
#pgls with averaged residuals or need to run brms 
#tapply to avg residuals per host species 

bmod=brm(richness~humanDisturbed+lcites+
           (1|obs)+(1|gr(label,cov=A)),
         chain=1,iter=3000,thin=1,warmup=1500,
         data=cdata$data,family=poisson,
         data2=list(A=cmatrix))

#brm example
#might want to look into if parasite random effects should be added
#for now just try pgls version with averaged residuals 
bmod=brm(jackknife~humanDisturbed+lcites+
           (1|obs)+(1|gr(label,cov=A)),
         chain=1,iter=3000,thin=1,warmup=1500,
         data=cdata$data,family=gaussian,
         data2=list(A=cmatrix))
