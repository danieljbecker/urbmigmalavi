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
library(plyr)

setwd("/home/verrett/R")

## load in MalAvi with traits
#setwd("~/Desktop/urbmigmalavi/flat files")
#setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
data=read.csv("MalAvi hosts with traits_cleaned_C_noLeuco.csv")

## load in citations
cdata=read.csv("urbmigmalavi host citations_noLeuco.csv")
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
data$UTI=as.numeric(data$UTI)

#fix migration
data$Migratory_status=revalue(data$Migration,
                     c("1"="not_full",
                       "2"="not_full",
                       "3"="full_migrant"))

## load BirdTree
#setwd("/Users/danielbecker/Desktop/urbmigmalavi/BirdTree")
#setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/BirdTree")
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
#ggplot(cdata$data[!is.na(cdata$data$urban),],
    #   aes(urban,richness))+geom_boxplot()+scale_y_log10()

#richness in human-disturbed
#ggplot(cdata$data[!is.na(cdata$data$humanDisturbed),],
     #  aes(humanDisturbed,richness))+geom_boxplot()+scale_y_log10()

## compute correlation matrix
cmatrix=vcv.phylo(cdata$phy,cor=T)

## disable
options(buildtools.check = function(action) TRUE)

## glmer (poisson with observation-level RE)
#mod=glmer(richness~urban*Migration+lcites+(1|obs),data=data,family=poisson)
#Anova(mod)

## glm.nb
#mod=glm.nb(richness~urban*Migration+lcites,data=data)
#Anova(mod)

## load in paco results
#setwd("~/Desktop/urbmigmalavi/flat files")
#setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
rdata=read.csv("PACo results_noLeuco.csv")
rdata$X=NULL
#collapsed_rdata=read.csv("collapsed PACo results_Passerellidae.csv")
#collapsed_rdata$X=NULL

## merge with traits
rdata=merge(rdata,data,by="tiplab",all.x=T)
#collapsed_rdata=merge(collapsed_rdata,data,by="tiplab",all.x=T)
library(dplyr)

#average jackknife for species and collapse species links
avg_jackknife <- tapply(rdata$jackknife, rdata$tiplab, mean)
collapsed_rdata = transform(rdata, average_jackknife = avg_jackknife[match(tiplab, names(avg_jackknife))])

#collapsed_rdata <- rdata %>%
#  group_by(tiplab) %>%
#  mutate(mean_jackknife = mean(jackknife)) %>%
#  distinct(tiplab, .keep_all = TRUE)

#cdata with paco data for pgls
#cdata=comparative.data(phy=htree,data=collapsed_rdata,names.col=tiplab,vcv=T,na.omit=F,warn.dropped=T)

## vis jacknife
#ggplot(collapsed_rdata[!is.na(collapsed_rdata$urban),],
       #aes(urban,average_jackknife))+geom_boxplot()

#ggplot(collapsed_rdata[!is.na(collapsed_rdata$humanDisturbed),],
       #aes(humanDisturbed,average_jackknife))+geom_boxplot()

#ggplot(collapsed_rdata[!is.na(collapsed_rdata$UTI),],
      # aes(UTI,average_jackknife))+geom_point()

#brm example
#might want to look into if parasite random effects should be added
#for now just try pgls version with averaged residuals 
bmod=brm(mean_jackknife~urban+Migratory_status+lcites+
        (1|gr(label,cov=A)),
         chains=4,iter=20000,thin=50,warmup=10000,
        control=list(adapt_delta = 0.98, max_treedepth=12),
         data=collapsed_rdata,family=gaussian,
         data2=list(A=cmatrix))
#this is currently returning a lot of errors (divergent transitions, exceeded max tree depth,
#low Bayesian Fraction of Missing Information, low bulk and tail ESS)

save.image(file='urbmigmod1_noLeuco.RData')

#hs=data.frame(bayestestR::hdi(bmod,ci=0.95))
#summary(bmod)
#mcmc_plot(bmod,type="trace") 