### This script is to calculate the cross-correlation between the temporal trends
### of subtidal predators identified to feed on the intertidal species and each of their
### subtidal preys. These were identified using binary connectance web with the GLOBI
### interaction data base from the scripts dataPrep_NOAA-intertidal_Web.R and
### NOAA-intertidal_interaction_web.R. Temporal trends in the intertidal are fitted
### in time-sries_JB.R. Temporal trends in the subtidal ar fitted in fit_NOAA_UV_TS_linkers.R.
### Contributor: Julien Beaulieu

# load packages
library(ggtext)
library(tidyr)
library(brms)
library(dplyr)
source("Scripts/pp_extract_function.R")
source("Scripts/function_genusSp_2_GSp.R")
library(ggplot2)
library(stringr)
library(ggpubr)

###________________________________________________________________
### make load intertidal model output, extract pp for each site and 
### put the averages in one DF ####################################
###________________________________________________________________


filenames <- list.files("time_series_outputs/new/final_UV_TS/")
mod_types <- c("I_zi", "I_hu", "I_hu", "I_zi","I_hu","I_zi","I_zi","I_hu",
              "I_zi", "I_zi", "I_zi")
sp_list_inter <- c("amphipoda","ascophyllum","Chondrus","cmaenas","fucus",
                   "Isopoda","littorina","mastocarpus","mytilus","nucella","semibalanus")

for (i in c(1:length(filenames))) {
  fit <- readRDS(paste("time_series_outputs/new/final_UV_TS/",filenames[i], sep = ""))
  
  pp <- extract(fit, model_type = mod_types[i])
  BC <- pp[[1]] %>% mutate(mean_pp = rowMeans(select(pp[[1]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "BC")
  NE <- pp[[2]] %>% mutate(mean_pp = rowMeans(select(pp[[2]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "NE")
  NW <- pp[[3]] %>% mutate(mean_pp = rowMeans(select(pp[[3]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "NW")
  SE <- pp[[4]] %>% mutate(mean_pp = rowMeans(select(pp[[4]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "SE")
  SW <- pp[[5]] %>% mutate(mean_pp = rowMeans(select(pp[[5]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "SW")
  zi_BC <- pp[[6]] %>% mutate(mean_pp = rowMeans(select(pp[[6]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "zi_BC")
  zi_NE <- pp[[7]] %>% mutate(mean_pp = rowMeans(select(pp[[7]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "zi_NE")
  zi_NW <- pp[[8]] %>% mutate(mean_pp = rowMeans(select(pp[[8]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "zi_NW")
  zi_SE <- pp[[9]] %>% mutate(mean_pp = rowMeans(select(pp[[9]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "zi_SE")
  zi_SW <- pp[[10]] %>% mutate(mean_pp = rowMeans(select(pp[[10]], -year))) %>% 
    select(c(year, mean_pp)) %>% mutate(site = "zi_SW")
  
  pp_all <- rbind(BC, NE, NW, SE, SW, zi_BC, zi_NE, zi_NW, zi_SE, zi_SW)
  
  assign(paste(sp_list_inter[i],"mean_pp", sep = "_"), pp_all)
  
  filenames[i]
  sp_list_inter[i]
}


### Load intraction data ###########################


sp_origin <- read.csv("data/species_origin_network.csv")
sub_sp <- sp_origin %>% filter(origin == "subidal") %>% mutate(specie = sub(" ","_",specie))
sub_sp <- as.vector(sub_sp$specie)
sp_inter <- sp_origin %>% filter(origin == "intertidal") %>% mutate(specie = sub(" ","_",specie))
sp_inter <- as.vector(sp_inter$specie)
interact <- read.csv("data/species_origin_directInter.csv") %>% 
  mutate(target = sub(" ","_",target)) %>% 
  # add missing important interaction based on Ellis et al 2007
  rbind(c("Cancer_borealis","eats","Nucella lapillus"),
        c("Cancer_borealis","eats","Littorina"),
        c("Cancer_borealis","eats","Carcinus maenas")) 
sp_list_sub <- interact %>% dplyr::select(target) %>% unique() %>% 
  filter(target %in% sub_sp) %>% mutate(target = gsub(" ","_",target))
sp_list_sub <- as.vector(sp_list_sub$target)


# add missing interactions

missing_inter <- read.csv("data/(Rest of) Subtidal Food Web - Interaction Data.csv") %>% 
  dplyr::select(-c(Paper.ID,Reference.DOI,X,localityName,
                   observationDateTime, Reference.ISBN, Notes, decimalLatitude,
                   decimalLongitude)) %>% 
  mutate(targetTaxonName = sub(" ", "_", targetTaxonName),
         sourceTaxonName = sub(" ", "_", sourceTaxonName)) %>% 
  #select trophic inter
  filter(interactionTypeName %in% c("preys on","is eaten by","preyed upon by","eats")) %>% 
  #convert everything in x eat y
  mutate(interactionTypeName = ifelse(interactionTypeName %in% c("eats","preys on"),
                                      "eats", "is eaten by")) %>%
  # here I am switching source and target because GLOBI data in form target eats source
  mutate(target = ifelse(interactionTypeName == "eats", sourceTaxonName, targetTaxonName),
         source = ifelse(interactionTypeName == "eats", targetTaxonName, sourceTaxonName),
         interaction_type = "eats") %>% 
  dplyr::select(target ,interaction_type, source) %>%
  # select only interaction with sub sp feeding ont int sp.
  filter(target %in% sub_sp,
         source %in% sp_inter) %>% 
  mutate(source = sub("_", " ", source)) %>% 
  unique()

interact <- unique(rbind(interact, missing_inter))

# remove useless stuff from the env.
rm(fit, NE, BC, NW,SE, SW, pp, pp_all, zi_NE, zi_BC, zi_NW, zi_SE, zi_SW)



###_____________________________________________________________________
### Make a loop to take each subtidal sp and calculate cross-correlation
### with preys that are in the key intertidal sp. ######################
###_____________________________________________________________________

# make cross_cor df to fill

cross_cor <- as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(cross_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")

for (i in c(1:length(sp_list_sub))) {
  
  # pull the subtidal specie model
  fit <- readRDS(paste("time_series_outputs/UV_TS_subtidal/",sp_list_sub[i],".RDS", sep = ""))
  
  # get the interactions data
  int_data <- interact %>% filter(target == sp_list_sub[i])
  
  
  ### check if any of the preys are Isopoda
  
  if(any(int_data$source == "Isopoda")){
    
    # get pp for the subtidal specie
    new_data <- Isopoda_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- Isopoda_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Isopoda_pp" = mean_pp)
    
    site_list <- unique(Isopoda_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Isopoda_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Isopoda") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  
  if(any(int_data$source == "Amphipoda")){
    
    # get pp for the subtidal specie
    new_data <- amphipoda_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- amphipoda_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Amphipoda_pp" = mean_pp)
    
    site_list <- unique(amphipoda_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Amphipoda_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Amphipoda") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  
  
  if(any(int_data$source == "Carcinus maenas")){
    
    # get pp for the subtidal specie
    new_data <- cmaenas_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- cmaenas_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Carcinus_maenas_pp" = mean_pp)
    
    site_list <- unique(cmaenas_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Carcinus_maenas_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Carcinus_maenas") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  
  if(any(int_data$source == "Littorina")){
    
    # get pp for the subtidal specie
    new_data <- littorina_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- littorina_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Littorina_pp" = mean_pp)
    
    site_list <- unique(littorina_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Littorina_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Littorina") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  if(any(int_data$source == "Mytilus edulis")){
    
    # get pp for the subtidal specie
    new_data <- mytilus_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- mytilus_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Mytilus_edulis_pp" = mean_pp)
    
    site_list <- unique(mytilus_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Mytilus_edulis_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Mytilus_edulis") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  
  
  
  if(any(int_data$source == "Fucus")){
    
    # get pp for the subtidal specie
    new_data <- fucus_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- fucus_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Fucus_pp" = mean_pp)
    
    site_list <- unique(fucus_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Fucus_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Fucus") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  if(any(int_data$source == "Chondrus crispus")){
    
    # get pp for the subtidal specie
    new_data <- Chondrus_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- Chondrus_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Chondrus_crispus_pp" = mean_pp)
    
    site_list <- unique(Chondrus_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Chondrus_crispus_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Chondrus_crispus") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  
  if(any(int_data$source == "Mastocarpus stellatus")){
    
    # get pp for the subtidal specie
    new_data <- mastocarpus_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- mastocarpus_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Mastocarpus_stellatus_pp" = mean_pp)
    
    site_list <- unique(mastocarpus_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Mastocarpus_stellatus_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Mastocarpus_stellatus") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  if(any(int_data$source == "Ascophyllum nodosum")){
    
    # get pp for the subtidal specie
    new_data <- ascophyllum_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- ascophyllum_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Ascophyllum_nodosum_pp" = mean_pp)
    
    site_list <- unique(ascophyllum_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Ascophyllum_nodosum_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Ascophyllum_nodosum") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  

  if(any(int_data$source == "Nucella lapillus")){
    
    # get pp for the subtidal specie
    new_data <- nucella_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- nucella_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Nucella_lapillus_pp" = mean_pp)
    
    site_list <- unique(nucella_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Nucella_lapillus_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Nucella_lapillus") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  
  
  
  if(any(int_data$source == "Semibalanus balanoides")){
    
    # get pp for the subtidal specie
    new_data <- semibalanus_mean_pp %>% rename("YEAR" = year)
    pred_data <- as.data.frame(predict(fit, newdata = new_data))
    dat_cor <- semibalanus_mean_pp %>% mutate(sub_pp = pred_data$Estimate) %>% 
      rename("Semibalanus_balanoides_pp" = mean_pp)
    
    site_list <- unique(semibalanus_mean_pp$site)
    
    # calculate cross correlation for each site
    for (j in c(1:length(site_list))) {
      
      
      temp <- dat_cor %>% filter(site == site_list[j])
      cc <- ccf(temp$Semibalanus_balanoides_pp, temp$sub_pp, plot = FALSE)
      
      temp_cor <- as.data.frame(matrix(nrow = length(cc[["acf"]]), ncol = 5))
      colnames(temp_cor) <- c("sub_sp", "inter_sp", "inter_site", "lag", "acf")
      
      temp_cor <- temp_cor %>% mutate(lag = cc[["lag"]]) %>% 
        mutate(correlation = cc[["acf"]]) %>% 
        mutate(sub_sp = sp_list_sub[i]) %>% 
        mutate(inter_sp = "Semibalanus_balanoides") %>% 
        mutate(site = site_list[j])
      
      cross_cor <- rbind(temp_cor, cross_cor)
    }
  }
  i
}

# re-formater pcq j'etais dans la lune

cross_cor <- cross_cor %>% select(-c(inter_site, acf))

# save table

write.csv(cross_cor, "time_series_outputs/cross_correlation_values_inter-sub2.csv",
          row.names = FALSE)

cross_cor <- read.csv("time_series_outputs/cross_correlation_values_inter-sub2.csv")

####___________________________________________________________________________
#### Make a heatmap with species and correlations #############################
####___________________________________________________________________________


### take only lag with max correlation

max_cor_lag <- cross_cor %>% unite("ID", c(sub_sp, inter_sp, site), remove = FALSE) %>% 
  group_by(ID) %>%
  mutate(abs_max_cor = max(abs(correlation))) %>% 
  ungroup() %>%
  filter(abs(correlation) == abs_max_cor) %>% 
  rename("Correlation" = correlation) %>% 
  select(-c(ID,lag,abs_max_cor)) %>% 
  unique() %>% unite("inter_sp_site", c(inter_sp, site), sep = "-") %>% 
  mutate(sub_sp = sub("_"," ", sub_sp)) %>% 
  separate(inter_sp_site, c("inter_sp","site"), sep = "-") %>% 
  mutate(inter_sp = sub("_"," ", inter_sp)) %>% 
  unite("inter_sp_site", c(inter_sp, site), sep = "-") 

### add a column to ID subtidal species that are impacted by at least one RS
### They were identified by looking at the species specific temporal trends


sub_sp_imp <- c(
  "Triglops murrayi",
  #"Scomber scombrus",#positively cross-correlated
  #"Pollachius virens",#positively cross-correlated
  "Peprilus triacanthus",
  "Prionotus carolinus",
  "Pandalus borealis",
  "Pseudopleuronectes americanus",
  #"Majidae",#positively cross-correlated
  #"Macrouridae",#positively cross-correlated
  "Lebbeus polaris",
  "Homarus americanus",
  #"Gadus morhua", #positively cross-correlated
  "Glyptocephalus cynoglossus",
  #"Cryptacanthodes maculatus",#positively cross-correlated
  #"Cyclopterus lumpus",#positively cross-correlated
  "Cancer borealis",
  #"Brosme brosme",#positively cross-correlated
  #"Amblyraja radiata",#positively cross-correlated
  "Anarhichas lupus"
)


sub_sp_imp <- unique(c(sub_sp_imp, "Alosa aestivalis",
                "Alosa pseudoharengus",
                "Aspidophoroides monopterygius",
                "Cancer borealis",
                "Chionoecetes opilio",
                "Chlorophthalmus agassizi",
                "Citharichthys arctifrons",
                "Clupea harengus",
                "Enchelyopus cimbrius",
                "Glyptocephalus cynoglossus",
                "Helicolenus dactylopterus",
                "Hippoglossoides platessoides",
                "Homarus americanus",
                "Melanogrammus aeglefinus",
                "Merluccius bilinearis",
                "Paralichthys oblongus",
                #"Pollachius virens", #positively cross-correlated
                "Pseudopleuronectes americanus",
                "Sebastes fasciatus",
                "Tautogolabrus adspersus",
                "Urophycis chuss",
                "Urophycis regia",
                "Urophycis tenuis"))


max_cor_lag <- max_cor_lag %>%
  mutate(sub_imp = ifelse(sub_sp %in% sub_sp_imp, "1", "0")) %>%
  genusSp2gSp(column = "sub_sp") %>% 
  mutate(sub_sp_label = ifelse(sub_imp == 1,
                               paste0("**", sub_sp, "**"),  # bold important species
                               sub_sp))

### plot heatmap


plot <- ggplot(data = max_cor_lag, aes(x = inter_sp_site, y = sub_sp_label, fill = Correlation)) +
  geom_tile(color = "black", lwd = 0.1, linetype = 1)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  ylab("Subtidal specie")+
  xlab("Intertidal specie - site")+
  theme_bw()+
  theme(axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        #axis.text.y = ggtext::element_markdown(size = 12),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 8),
        panel.grid = element_line(color = "black",linewidth = 0.5))



plot <- ggplot(
  data = max_cor_lag,
  aes(x = inter_sp_site, y = sub_sp, fill = Correlation)
) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_y_discrete(
    labels = setNames(max_cor_lag$sub_sp_label, max_cor_lag$sub_sp)
  ) +
  ylab("Subtidal species") +
  xlab("Intertidal species - site") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_markdown(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 8),
    panel.grid = element_line(color = "black", linewidth = 0.5)
  )

plot

ggsave("time_series_outputs/figures/cross_corr_heatmap.tiff",
       plot = plot, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")



### make another figure separating the zi trend in another panel


max_cor_zi <- max_cor_lag %>% 
  filter(str_detect(inter_sp_site, "zi"))

plot_zi <- ggplot(data = max_cor_zi, aes(x = inter_sp_site, y = sub_sp, fill = Correlation)) +
  geom_tile(color = "black", lwd = 0.1, linetype = 1)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  coord_fixed()+
  ylab("Subtidal specie")+
  xlab("Intertidal specie - site")+
  theme_bw()+
  theme(axis.text=element_text(size=9),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 7),
        panel.grid = element_line(color = "black",linewidth = 0.5))

plot_zi



max_cor_mu <- max_cor_lag %>% 
  filter(str_detect(inter_sp_site, "zi", negate = TRUE)) %>%
  genusSp2gSp(column = "sub_sp")
  

plot_mu <- ggplot(data = max_cor_mu, aes(x = inter_sp_site, y = sub_sp, fill = Correlation)) +
  geom_tile(color = "black", lwd = 0.1, linetype = 1)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  scale_y_discrete(labels = setNames(max_cor_mu$sub_sp_label, max_cor_mu$sub_sp)) +
  ylab("Subtidal species")+
  xlab("Intertidal species - sites")+
  theme_bw()+

  theme(axis.text=element_text(size=8.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 10),
        axis.text.y = element_markdown(size = 8.5),
        legend.title = element_text(size = 12),
        panel.grid = element_line(color = "black",linewidth = 0.1))#,
        #legend.position = "none")

plot_mu

ggsave("time_series_outputs/figures/cross_corr_heatmap_mu2.tiff",
       plot = plot_mu, device = "tiff", scale = 2, 
       width = 1500, height = 1200, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
## make multipanel figure


fig <- ggarrange(plot_mu, plot_zi,
          labels = c("a) Count model", "b) Probability model (zi)"),
          hjust = -0.2, vjust = -0.075,
          font.label = list(size = 15, face = "bold.italic"),
          nrow = 1, ncol = 2) %>% 
  annotate_figure(top = text_grob("", color = "white", size = 20), 
  )
fig

ggsave("time_series_outputs/figures/cross_corr_heatmap.tiff",
       plot = fig, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")


### make a multipanel figure with the count model and the network
### from the NOAA-intertidal_interaction_web.R


fig <- ggarrange(plot2, plot_mu,
                 labels = c("a) Binary connectance web", "b) Cross-correlation"),
                 hjust = -0.2, vjust = -0.075,
                 font.label = list(size = 15, face = "bold.italic"),
                 nrow = 1, ncol = 2, widths = c(0.75,1)) %>% 
  annotate_figure(top = text_grob("", color = "white", size = 20), 
  )
fig

ggsave("time_series_outputs/figures/web_cross-corr_option2.tiff",
       plot = fig, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

