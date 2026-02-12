# This script is to check and plot multivariate time series

#load packages

library(ggplot2)
library(tidyverse)
library(rstan)
library(devtools)
library(brms)
library(rstan)
library(rstantools)
library(mgcv)
library(schoenberg)
library(coda)
library(bayesplot)
library(loo)

# load data
dat <- read_csv("data/full_dataset_for_SEMs.csv", guess_max = 3000) #data produced by SEM_dataset_assembly.R

# Add rows for sum littorina and fucus 

dat$sum_littorina <- rowSums(dat[,c("littorina_littorea", "littorina_obtusata", "littorina_saxatilis")])
dat$sum_fucus <- rowSums(dat[,c("fucus_vesiculosus","fucus_distichus")])
dat$intertidal_transect <- as.factor(dat$intertidal_transect)

# load models results
noaa <- readRDS("time_series_outputs/noaa_pcoa_GS_3.rds")


GS <- readRDS("time_series_outputs/output/count_pcoa_GS_10.rds")
GS_count <- add_criterion(GS, c("loo","bayes_R2"))
GI <- readRDS("time_series_outputs/output/count_pcoa_GI_10.rds")
GI_count <- add_criterion(GI, c("loo","bayes_R2"))
I <- readRDS("time_series_outputs/output/count_pcoa_I_10.rds")
I_count <- add_criterion(I,c("loo","bayes_R2"))
count_4 <- readRDS("time_series_outputs/new/count_pcoa_mix_9.rds")

GS <- readRDS("time_series_outputs/output/cover_pcoa_GS_6.rds") 
GI <- readRDS("time_series_outputs/output/cover_pcoa_GI_10.rds")## Best model
I <- readRDS("time_series_outputs/output/cover_pcoa_I_10.rds")



# check convergence

mcmc_plot(GS, type = "trace")
mcmc_plot(GI, type = "trace")
mcmc_plot(I, type = "trace")
mcmc_plot(count_4, type = "trace")

# check the general fit

pp_check(GS, resp = "MDS1")
pp_check(GS, resp = "MDS2")
pp_check(GS, resp = "MDS3")

pp_check(GI, resp = "MDS1")
pp_check(GI, resp = "MDS2")
pp_check(GI, resp = "MDS3")

pp_check(I, resp = "MDS1")
pp_check(I, resp = "MDS2")
pp_check(I, resp = "MDS3")

pp_check(count_4, resp = "MDS1")
pp_check(count_4, resp = "MDS2")
pp_check(count_4, resp = "MDS3")

# check the fit per site

pp_check(GS, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS1") 
pp_check(GI, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS1")
pp_check(I, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS1") 

pp_check(GS, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS2") 
pp_check(GI, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS2")
pp_check(I, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS2")

pp_check(GS, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS3") 
pp_check(GI, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS3")
pp_check(I, type = "ribbon_grouped", group = "SITE", nsamples = 25, resp = "MDS3")


# cross validation

GS <- add_criterion(GS, "loo_subsample")
GI <- add_criterion(GI, "loo_subsample")
I <- add_criterion(I, "loo_subsample") 

GS$criteria$loo
GS$criteria$waic
GS$criteria$bayes_R2

I$criteria$loo
I$criteria$waic
I$criteria$bayes_R2

# visualise pp
msms_noaa <- conditional_smooths(noaa)
plot(msms_noaa)

msms_GS <- conditional_smooths(GS)
plot(msms_GS)

msms_GI <- conditional_smooths(GI)
plot(msms_GI)

msms_I <- conditional_smooths(I)
plot(msms_I)

#Posterior predictions
p_noaa <- conditional_effects(noaa)
p_noaa

p_GS <- conditional_effects(GS)
p_GS

p_GI <- conditional_effects(GI)
p_GI


p_I <- conditional_effects(I) 
p_I
