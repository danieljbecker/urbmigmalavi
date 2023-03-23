## MalAvi, urbanization, migration
## 01_MalAvi clean
## danbeck@ou.edu
## last updated 3/23/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## libraries
library(malaviR)
library(plyr)

## load MalAvi host-lineage data
data=extract_table("Hosts and Sites Table")

## trim columns
data=data[c("Lineage_Name","Alt_Lineage_Name","parasiteGenus","species")]

## make unique
data$id=paste(data$Lineage_Name,data$species)

## remove dups
data=data[!duplicated(data$id),]
data$id=NULL

## flag hybrids
x=strsplit(data$species," ")
data$length=sapply(x,length)

## flag genus only
x=sapply(x,function(x) x[2])
data$sp=x

## remove hybrids and genus only
data=data[!data$length>2,]
data=data[!data$sp%in%c("sp","sp.","spp","spp."),]

## clean
rm(x)
data$length=NULL
data$sp=NULL

## save link column
data$mal_original=data$species

## match taxonomy
tax=taxonomy

## merge in with data
data=merge(data,tax,by="species",all.x=T)

## fix Jetz
data$Jetz.species=as.character(data$Jetz.species)
data$Jetz.species=ifelse(data$Jetz.species=="",NA,data$Jetz.species)

## if Jetz match, use that name
data$tip=ifelse(is.na(data$Jetz.species),data$species,data$Jetz.species)

## fix tip
data$tip=gsub("_"," ",data$tip)

## load BirdTree taxonomy for backbone
setwd("~/OneDrive - University of Oklahoma/Becker Lab/Datasets/BirdTree")
pnames=read.csv("BLIOCPhyloMasterTax.csv")

## check mismatch
mis=setdiff(data$tip,pnames$Scientific) ## 80 mismatch

## manually edit remaining 80 mismatch
data$tip=revalue(data$tip,
                 c("Acritillas indica"="Iole indica",
                   "Actitis macularia"="Actitis macularius",
                   "Anisognathus flavinuchus"="Anisognathus somptuosus",
                   "Anser canagica"="Chen canagica",
                   "Anthropoides paradiseus"="Grus paradisea",
                   "Anthus cinnamomeus"="Anthus novaeseelandiae",
                   "Bubo scandiacus"="Bubo scandiaca",
                   "Cantorchilus longirostris"="Thryothorus longirostris",
                   "Cantorchilus superciliaris"="Thryothorus superciliaris",
                   "Cardellina canadensis"="Wilsonia canadensis",
                   "Chelidorhynx hypoxanthus"="Chelidorhynx hypoxantha",
                   "Chloris sinica"="Carduelis sinica",
                   "Chlorospingus flavopectus"="Chlorospingus ophthalmicus",
                   "Chroicocephalus ridibundus"="Larus ridibundus",
                   "Cichlherminia lherminieri"="Turdus lherminieri",
                   "Clanga pomarina"="Aquila pomarina",
                   "Corythornis madagascariensis"="Ceyx madagascariensis",
                   "Criniger pallidus"="Alophoixus pallidus",
                   "Cyanoloxia brissonii"="Cyanocompsa brissonii",
                   "Cyornis pallidipes"="Cyornis pallipes",
                   "Cypseloides rutilus"="Streptoprocne rutila",
                   "Daptrius americanus"="Ibycter americanus",
                   "Diglossopis glauca"="Diglossa glauca",
                   "Doryfera ludoviciae"="Doryfera ludovicae",
                   "Elaenia chilensis"="Elaenia albiceps",
                   "Erythrogenys erythrogenys"="Pomatorhinus erythrogenys",
                   "Euaegotheles insignis"="Aegotheles insignis",
                   "Falcipennis canadensis"="Dendragapus canadensis",
                   "Geokichla citrina"="Zoothera citrina",
                   "Geospizopsis unicolor"="Phrygilus unicolor",
                   "Griseotyrannus aurantioatrocristatus"="Empidonomus aurantioatrocristatus",
                   "Gymnoris pyrgita"="Petronia pyrgita",
                   "Hartlaubius auratus"="Saroglossa aurata",
                   "Henicorhina anachoreta"="Henicorhina leucophrys",
                   "Hesperiphona vespertina"="Coccothraustes vespertinus",
                   "Ispidina picta"="Ceyx pictus",
                   "Ketupa ketupa"="Ketupa ketupu",
                   "Kittacincla malabarica"="Copsychus malabaricus",
                   "Machlolophus spilonotus"="Parus spilonotus",
                   "Malurus assimilis"="Malurus lamberti",
                   "Melanitta americana"="Melanitta nigra",
                   "Melanitta stejnegeri"="Melanitta fusca",
                   "Microscelis amaurotis"="Ixos amaurotis",
                   "Mixornis gularis"="Macronous gularis",
                   "Momotus lessonii"="Momotus aequatorialis",
                   "Montecincla jerdoni"="Strophocincla cachinnans",
                   "Montecincla meridionalis"="Strophocincla fairbanki",
                   "Myiothlypis flaveola"="Basileuterus flaveolus",
                   "Nyctipolus nigrescens"="Caprimulgus nigrescens",
                   "Oedistoma iliolophus"="Toxorhamphus iliolophus",
                   "Origma murina"="Crateroscelis murina",
                   "Origma robusta"="Crateroscelis robusta",
                   "Parus minor"="Parus major",
                   "Peneothello sigillata"="Peneothello sigillatus",
                   "Pheugopedius eisenmanni"="Thryothorus eisenmanni",
                   "Philohydor lictor"="Pitangus lictor",
                   "Piculus rivolii"="Colaptes rivolii",
                   "Piculus rubiginosus"="Colaptes rubiginosus",
                   "Poecile sclateri"="Parus sclateri",
                   "Psittiparus gularis"="Paradoxornis gularis",
                   "Rhynchospiza stolzmanni"="Aimophila stolzmanni",
                   "Rupornis magnirostris"="Buteo magnirostris",
                   "Sciaphylax hemimelaena"="Myrmeciza hemimelaena",
                   "Seicercus xanthodryas"="Phylloscopus borealis",
                   "Setophaga citrina"="Wilsonia citrina",
                   "Sholicola ashambuensis"="Myiomela albiventris",
                   "Sholicola major"="Myiomela major",
                   "Silvicultrix frontalis"="Ochthoeca frontalis",
                   "Silvicultrix jelskii"="Ochthoeca jelskii",
                   "Spizaetus cirrhatus"="Nisaetus cirrhatus",
                   "Sporophila angolensis"="Oryzoborus angolensis",
                   "Stomiopera unicolor"="Lichenostomus unicolor",
                   "Syndactyla ucayalae"="Simoxenops ucayalae",
                   "Thamnophilus bernardi"="Sakesphorus bernardi",
                   "Tringa brevipes"="Heteroscelus brevipes",
                   "Trochalopteron cachinnans"="Strophocincla cachinnans",
                   "Trochalopteron fairbanki"="Strophocincla fairbanki",
                   "Turdus eunomus"="Turdus naumanni",
                   "Tyto furcata"="Tyto alba",
                   "Zoothera aurea"="Zoothera dauma"))

## recheck
mis=setdiff(data$tip,pnames$Scientific)
rm(mis)

## merge
pnames$tip=pnames$Scientific
pnames$Scientific=NULL
data=merge(data,pnames,by="tip",all.x=T)
rm(pnames)

## make id
data$id=paste(data$Lineage_Name,data$tip)

## remove dups
dups=names(which(table(data$id)>1))
tmp=data[data$id%in%dups,]

## remove from data
data=data[!data$id%in%dups,]

## preferentially sort
lset=list()
for(i in 1:length(unique(tmp$id))){
  
  ## subset
  set=tmp[tmp$id==unique(tmp$id)[i],]
  
  ## sort
  set=set[order(set$Alt_Lineage_Name),]
  set=set[!duplicated(set$id),]
  lset[[i]]=set
}

## combine
tmp=do.call(rbind.data.frame,lset)

## remerge with data
data=rbind.data.frame(data,tmp)
rm(set,tmp,dups,i,lset,tax)

## save edgelist
setwd("~/Desktop/urbmigmalavi/flat files")
edge=data
write.csv(edge,"MalAvi edgelist.csv")