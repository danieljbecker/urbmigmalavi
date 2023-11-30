library(stringr)
library(dplyr)

setwd("/Users/taylorverrett/Documents/GitHub/urbmigmalavi/flat files")
edgelist=read.csv("Malavi edgelist.csv")

split <- str_split_fixed(edgelist$species, " ", 2)
split = as.data.frame(split)
edgelist$genus = split$V1
rm(split)

sparrow_edgelist <- edgelist %>%
  filter(genus %in% c("Oreothraupis", "Chlorospingus", "Rhynchospiza", "Peucaea",
                      "Ammodramus", "Arremonops", "Amphispizopsis", "Amphispiza",
                      "Chondestes", "Calamospiza", "Spizella", "Arremon",
                      "Passerella", "Spizelloides", "Junco", "Zonotrichia",
                      "Artemisiospiza", "Oriturus", "Pooecetes", "Ammospiza",
                      "Centronyx", "Passerculus", "Xenospiza", "Melospiza",
                      "Pezopetes", "Torreornis", "Melozone", "Aimophila",
                      "Pipilo", "Atlapetes"))

write.csv(sparrow_edgelist,"sparrows only MalAvi edgelist.csv")