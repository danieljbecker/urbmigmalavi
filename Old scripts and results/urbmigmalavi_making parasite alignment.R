#install.packages("devtools")
#install.packages("rlang")
#devtools::install_github("vincenzoaellis/malaviR", build_vignettes = TRUE, force=TRUE)

library(malaviR)

#sequence table
data=extract_table("Grand Lineage Summary")

#sequences and name only
data=data[c("Lineage_Name","sequence")]

#dan's file
cleaned_lineages=read.csv("Malavi lineages cleaned.csv")

#delete extras
mis=setdiff(data$Lineage_Name,cleaned_lineages$Lineage_Name)
cleaned_sequences=data[!data$Lineage_Name%in%mis,]

#no idea why there is 1 extra observation, it insists the lineage names match
mis=setdiff(data$Lineage_Name,cleaned_lineages$Lineage_Name)

#save csv
write.csv(cleaned_sequences, file="/Users/taylorverrett/Desktop/Projects/cleaned sequences.csv", row.names=FALSE)

#convert to fasta
library (tidyverse)
library(devtools)

source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

TabularToFasta("cleaned sequences.csv")

