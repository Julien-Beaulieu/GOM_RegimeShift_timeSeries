# This script is to compare the species present in the NOAA bottom-trawl surveys 
# to those observed in the Shouls subtidal survey (keen data)
# Author: Julien Beaulieu

# load packages
library(dplyr)
source("Scripts/format_sp_name_for_globi.R")
library(tidyr)

noaa_sp <- read.csv("time_series_outputs/ordination_output/score_noaa_specie.csv") %>% 
  select(species) %>% 
  fix_sp_name() %>% 
  mutate(specie = ifelse(specie == "Sum cancer", "Cancer borealis", 
                         ifelse(specie == "Sum urophycis", "Urophycis regia",
                         specie)))
  

keen_sp <- read.csv("data/KEEN_data.csv") %>% 
  select(SPECIES) %>% 
  unique() %>% na.omit()


## identify species overlapping and not overlapping

which(noaa_sp$specie %in% keen_sp$SPECIES) ## ISHHHH!!!!!

length(which(noaa_sp$specie %in% keen_sp$SPECIES))/nrow(noaa_sp)*100

noaa_sp$specie[which(noaa_sp$specie %in% keen_sp$SPECIES)]
