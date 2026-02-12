### This script is to fit a GAM on each subtidal specie from the NOAA
### bottom trawl survey that has a direct interaction with intertidal taxa.
### These species are identified using binary connectance web with the GLOBI
### interaction data base from the scripts dataPrep_NOAA-intertidal_Web.R and
### NOAA-intertidal_interaction_web.R
# contributor: Julien Beaulieu

## load packages
library(tidyverse)
library(devtools)
library(brms)
library(rstan)
library(rstantools)
library(mgcv)
library(schoenberg)
library(coda)
library(bayesplot)
library(janitor)
library(reshape2)
library(vegan)
library(data.table)
library(ggrepel)
source("Scripts/format_sp_name_for_globi.R")

### Load and prep data

noaa <- read_csv("data/NOAA_trawl_data.csv") %>%
  filter(scaled_med_individ_PUE != "NaN") %>% #there are some NaNs in the data but don't include any of our species of interest
  filter(scaled_med_tot_mass_PUE != "NaN") %>%
  mutate(scaled_med_individ_PUE = as.numeric(scaled_med_individ_PUE)) %>%
  mutate(scaled_med_tot_mass_PUE = as.numeric(scaled_med_tot_mass_PUE)) %>% 
  mutate(sciname = make_clean_names(SCINAME)) %>% fix_sp_name(column = 10) %>% 
  mutate(specie = gsub(" ","_",specie))

#-Convert to wide dataframe-#
noaa_wide <- dcast(noaa, YEAR ~ specie, value.var = "scaled_med_tot_mass_PUE", mean)

#-NaNs are years where no individuals of that species was detected, will replace with 0-#

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
noaa_wide[is.nan(noaa_wide)] <- 0

# sum species with the same genus
noaa_wide$sum_alosa <- rowSums(noaa_wide[,c("Alosa_aestivalis", "Alosa_mediocris", "Alosa_pseudoharengus","Alosa_sapidissima")])
noaa_wide$sum_ammodytes <- rowSums(noaa_wide[,c("Ammodytes_americanus","Ammodytes_dubius")])
noaa_wide$sum_argentina <- rowSums(noaa_wide[,c("Argentina_silus","Argentina_striata")])
noaa_wide$sum_cancer <- rowSums(noaa_wide[,c("Cancer_borealis","Cancer_irroratus")])
noaa_wide$sum_urophycis <- rowSums(noaa_wide[c("Urophycis_chesteri", "Urophycis_chuss","Urophycis_regia",
                                               "Urophycis","Urophycis_tenuis")])


sp_origin <- read.csv("data/species_origin_network.csv")
sub_sp <- sp_origin %>% filter(origin == "subidal")
sub_sp <- as.vector(sub_sp$specie)
interact <- read.csv("data/species_origin_directInter.csv")
sp_list <- interact %>% dplyr::select(target) %>% unique() %>% 
  filter(target %in% sub_sp) %>% mutate(target = gsub(" ","_",target))
sp_list <- as.vector(sp_list$target)

#### Make a loop that fit and saves TS for species in sp_list

for (i in c(1:length(sp_list))) {
  temp <- noaa_wide %>% dplyr::select(YEAR,sp_list[i]) %>% 
    rename("specie" = sp_list[i])
  
  mod <- bf(specie ~ s(YEAR), family = gaussian)
  get_prior(mod, temp)
  
  prior <- c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0.7, 2.5), class = "Intercept", lb = 0),
    prior(student_t(3, 0, 2.5), class = "sds", lb = 0),
    prior(student_t(3, 0, 2.5), class = "sigma", lb = 0)
    #prior(gamma(0.01,0.01), class = "shape", lb = 0)
  )

  fit <- brm(
    formula = mod,
    prior = prior,
    data = temp,
    cores = 4,
    chains = 4,
    iter = 10000,
    warmup = 5000,
    control = list(adapt_delta = 0.99)
  )
  
  fit <- add_criterion(fit, c("loo","bayes_R2","loo_R2"))
  
   # mcmc_plot(fit, type = "trace")
   # mcmc_plot(fit, type = "hist")
   # pp_check(fit)

  saveRDS(fit, file = paste("time_series_outputs/UV_TS_subtidal/",sp_list[i],".RDS", sep = ""))
}
 
