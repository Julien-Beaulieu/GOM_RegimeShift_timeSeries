## This script is to pull the interactions of the NOAA subtidal sp. 
## with the intertidal species to identify key linkers. The binary connectance web
## using this data is made in the script NOAA-intertidal_interaction_web.R.
## This is the first step to propos a mechanism for RS propagation.
## Contributor: Julien Beaulieu

### Load the data
NOAA <- read.csv( "time_series_outputs/ordination_output/score_noaa_specie.csv")
NOAA <- NOAA[,-1]
colnames(NOAA)[1] <- "specie"
cover <- read.csv("time_series_outputs/ordination_output/spp_score_cover.csv")
cover <- cover[,-1]
count <- read.csv("time_series_outputs/ordination_output/spp_score_count.csv")
count <- count[,-1]

### Load packages

library(igraph)
library(dplyr)
library(tidyverse)
library(ggnetwork)
library(ggplot2)
library(NetIndices)
library(rglobi)

#load function to fix sp name to get globi working
source("Scripts/format_sp_name_for_globi.R")

### Prep intertidal sp. list ################

int_sp <- rbind(cover, count) %>% dplyr::select(species) %>%
  mutate(species = ifelse(species == "isopods","isopoda",
                          ifelse(species == "nematodes","nematoda",
                                 ifelse(species == "sum_littorina","littorina",
                                        ifelse(species == "corallinne", "corallinales",
                                               ifelse(species == "carcinus_maenus","carcinus_maenas",species)))))) %>% 
  fix_sp_name()


### Make subtidal sp list and pull interactions ################

mis_sp <- as.data.frame(c("alosa_aestivalis", "alosa_mediocris", "alosa_pseudoharengus","alosa_sapidissima",
            "ammodytes_americanus","ammodytes_dubius","argentina_silus","argentina_striata",
            "cancer_borealis","cancer_irroratus","urophycis_chesteri", "urophycis_chuss","urophycis_regia",
            "urophycis_sp","urophycis_tenuis","chanceon_quinquedens")) # chanceon = Geryon
colnames(mis_sp) <- "specie"

### prep sp list
sub_sp <- NOAA %>% mutate(specie = ifelse(specie == "sum_alosa", "Alosa_sp",
                          ifelse(specie == "sum_ammodytes", "Ammodytes_sp",
                                 ifelse(specie == "sum_argentina", "Argentina_sp",
                                        ifelse(specie == "sum_cancer", "Cancer_sp",
                                               ifelse(specie == "sum_urophycis", "Urophycis_sp",
                                                      ifelse(specie == "lycenchelys_verrilli","lycenchelys_verrillii",specie))))))) %>%
  dplyr::select(specie) %>%
  rbind(mis_sp) %>% 
  fix_sp_name() %>%  
  # remove sum of stuff and replace by species present
  filter(!(specie %in% c("Alosa","Ammodytes","Argentina","Cancer","Urophycis")))
  
### pull from GLOBI
sub_sp <- as.vector(sub_sp$specie)
int_sp <- as.vector(int_sp$specie)

all_taxa <- c(sub_sp,int_sp)

interact <- get_interactions_by_taxa(all_taxa[1])
interact$taxon_called = all_taxa[1]
n=0
for (i in c(2:length(all_taxa))) {
  temp <- get_interactions_by_taxa(all_taxa[i])
  if(nrow(temp) >= 1){temp$taxon_called = all_taxa[i]}
  if(nrow(temp) < 1){
    print(all_taxa[i])
    n = n+1
    print(n)}
  interact <- rbind(interact, temp)
}


### Add missing interactions

missing_inter <- read.csv("data/(Rest of) Subtidal Food Web - Interaction Data.csv") %>% 
  dplyr::select(-c(Paper.ID,Reference.DOI,X,localityName,
                   observationDateTime, Reference.ISBN, Notes, decimalLatitude,
                   decimalLongitude)) %>% 
  #select trophic inter
  filter(interactionTypeName %in% c("preys on","is eaten by","preyed upon by","eats")) %>% 
  #convert everything in x eat y
  mutate(interactionTypeName = ifelse(interactionTypeName %in% c("eats","preys on"),
                                      "eats", "is eaten by")) %>%
  # here I am switching source and target because GLOBI data in form target eats source
  mutate(target_taxon_name = ifelse(interactionTypeName == "eats", sourceTaxonName, targetTaxonName),
         source_taxon_name = ifelse(interactionTypeName == "eats", targetTaxonName, sourceTaxonName),
         interaction_type = "eats") %>% 
  dplyr::select(target_taxon_name ,interaction_type, source_taxon_name) %>% 
  mutate(source_taxon_external_id = NA,
         source_taxon_path = NA,
         source_specimen_life_stage = NA,
         target_taxon_external_id = NA,
         target_taxon_path = NA,
         target_specimen_life_stage = NA,
         latitude = NA,
         longitude = NA,
         study_citation = NA,
         study_source_citation = NA,
         taxon_called = source_taxon_name) %>% 
  unique()


interact <- rbind(interact, missing_inter)

### Select subtidal species feeding on any of the intertidal sp. #############

inter_cl <- interact %>% 
  dplyr::select(source_taxon_name,source_taxon_path,interaction_type,
         target_taxon_name,target_taxon_path, taxon_called) %>% 
  filter(interaction_type %in% c("eats","preysOn","kills")) %>% 
  filter(str_detect(target_taxon_name, paste(int_sp, collapse = "|"))|
           str_detect(target_taxon_path, paste(int_sp, collapse = "|")))

# identify which intertidal taxa is the target 
for (i in c(1:nrow(inter_cl))) {
  for (j in c(1:length(int_sp))) {
    if (str_detect(inter_cl$target_taxon_name[i], int_sp[j])) {
      inter_cl$target[i] = int_sp[j]
    }
    if (is.na(inter_cl$target_taxon_path[i]) == FALSE) {
      if (str_detect(inter_cl$target_taxon_path[i], int_sp[j])) {
        inter_cl$target[i] = int_sp[j]
      }
    }
  }
}


direct_inter <- inter_cl %>% dplyr::select(c(taxon_called,interaction_type,target)) %>% 
  rename("source" = target, "target" = taxon_called)%>% 
  mutate(interaction_type = "eats") %>% unique()



write.csv(direct_inter, file = "data/direct_interaction_int-sub.csv",
          row.names = FALSE)

### select species interacting with species interacting with intertidal species #########

indir_dir_sp <- c(int_sp, direct_inter$target)

indirect_direct_inter <- interact %>% 
  dplyr::select(source_taxon_name,source_taxon_path,interaction_type,
                target_taxon_name,target_taxon_path, taxon_called) %>% 
  filter(interaction_type %in% c("eats","preysOn","kills")) %>% 
  filter(str_detect(target_taxon_name, paste(indir_dir_sp, collapse = "|"))|
           str_detect(target_taxon_path, paste(indir_dir_sp, collapse = "|"))) 

# identify which intertidal taxa is the target 
for (i in c(1:nrow(indirect_direct_inter))) {
  for (j in c(1:length(indir_dir_sp))) {
    if (str_detect(indirect_direct_inter$target_taxon_name[i], indir_dir_sp[j])) {
      indirect_direct_inter$target[i] = indir_dir_sp[j]
    }
    if (is.na(indirect_direct_inter$target_taxon_path[i]) == FALSE) {
      if (str_detect(indirect_direct_inter$target_taxon_path[i], indir_dir_sp[j])) {
        indirect_direct_inter$target[i] = indir_dir_sp[j]
      }
    }
  }
}

indirect_direct_inter_cl <- indirect_direct_inter %>% 
  dplyr::select(c(taxon_called,interaction_type,target)) %>% 
  rename("source" = target, "target" = taxon_called)%>% 
  mutate(interaction_type = "eats") %>% unique()
  

write.csv(indirect_direct_inter_cl, file = "data/indirect_direct_interaction_int-sub.csv",
          row.names = FALSE)



### make a dataframe with all taxa and their origin (int sub) ################################

#verify there is no overlap

for (i in c(1:length(sub_sp))) {
  for (j in c(1:length(int_sp))) {
    if(sub_sp[i] == int_sp[j]){print(i,j,sub_sp[i])}
  }
}


int_sp <- as.data.frame(int_sp) %>% mutate(origin = "intertidal") %>% 
  rename("specie" = int_sp)

all_sp <- as.data.frame(sub_sp) %>% mutate(origin = "subidal") %>% 
  rename("specie" = sub_sp) %>% 
  rbind(int_sp)

write.csv(all_sp, "data/species_origin_network.csv", row.names = FALSE)

