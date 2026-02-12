# Time series
# This script is to do model selection on the principal organism from the SEMs
# Contributor: Julien Beaulieu

# 0: load packages and data  #################################################

library(ggplot2)
library(tidyverse)
library(rstan)
library(devtools)
library(brms)
library(rstan)
library(rstantools)
library(schoenberg)
library(coda)
library(bayesplot)
library(loo)

# load data
dat <- read_csv("data/full_dataset_for_SEMs.csv", guess_max = 3000) #data produced by SEM_dataset_assembly.R

# Add rows for sum littorina and fucus 

dat$sum_littorina <- rowSums(dat[,c("littorina_littorea", "littorina_obtusata","littorina_saxatilis")])
dat$sum_fucus <- rowSums(dat[,c("fucus_vesiculosus","fucus_distichus")])
dat$intertidal_transect <- as.factor(dat$intertidal_transect)


# 1: load models results #######################################################
#_________________________________________________________________________


# load models results
#mytilus_GI <- readRDS("time_series_outputs/mytilusmod_TS_GI.rds")

#GS <- readRDS("time_series_outputs/cover_pcoa_GS_4.rds")
#I <- readRDS("time_series_outputs/cover_pcoa_I_4.rds")

#GS <- readRDS("time_series_outputs/output/mytilus_GS_3.rds")
#GI <- readRDS("time_series_outputs/output/mytilus_GI_3.rds")
#I <- readRDS("time_series_outputs/output/mytilus_I_3.rds") #best model

#GS <- readRDS("time_series_outputs/output/semibalanus_GS_3.rds") #best model
#GI <- readRDS("time_series_outputs/output/semibalanus_GI_3.rds")
#I <- readRDS("time_series_outputs/output/semibalanus_I_3.rds")


#GS <- readRDS("time_series_outputs/output/nucella_GS_3.rds")
#GI <- readRDS("time_series_outputs/output/nucella_GI_3.rds")
#I <- readRDS("time_series_outputs/output/nucella_I_3.rds") #best model

#GS <- readRDS("time_series_outputs/output/mastocarpus_GS_3.rds")
#GI <- readRDS("time_series_outputs/output/mastocarpus_GI_3.rds")#best model
#I <- readRDS("time_series_outputs/output/mastocarpus_I_3.rds")

#GS <- readRDS("time_series_outputs/output/littorina_GS_3.rds")
#GI <- readRDS("time_series_outputs/output/littorina_GI_3.rds")
#I <- readRDS("time_series_outputs/output/littorina_I_3.rds")

#GS <- readRDS("time_series_outputs/output/fucus_GS_3.rds") #best model
#GI <- readRDS("time_series_outputs/output/fucus_GI_3.rds")
#I <- readRDS("time_series_outputs/output/fucus_I_3.rds")

#GS <- readRDS("time_series_outputs/output/cmenas_GS_3.rds")
#GI <- readRDS("time_series_outputs/output/cmenas_GI_3.rds")
#I <- readRDS("time_series_outputs/output/cmenas_I_3.rds") #best model


#GS <- readRDS("time_series_outputs/output/chondrus_GS_3.rds")
#GI <- readRDS("time_series_outputs/output/chondrus_GI_3.rds")
#I <- readRDS("time_series_outputs/output/chondrus_I_3.rds") # best model

#GS <- readRDS("time_series_outputs/output/ascophyllum_GS_3.rds") #best model
#GI <- readRDS("time_series_outputs/output/ascophyllum_GI_3.rds")
#I <- readRDS("time_series_outputs/output/ascophyllum_I_3.rds")

# 2nd difference with wave exposure
# 
# mytilus_2_WE <- readRDS("time_series_outputs/new/mytilus_I_8.rds")
# semibalanus_2_WE <- readRDS("time_series_outputs/new/semibalanus_GS_8.rds")
# nucella_2_WE <- readRDS("time_series_outputs/new/nucella_I_8.rds")
# littorina_2_WE <- readRDS("time_series_outputs/new/littorina_I_8.rds")
# mastocarpus_2_WE <- readRDS("time_series_outputs/new/mastocarpus_GI_8.rds")
# fucus_2_WE <- readRDS("time_series_outputs/new/fucus_GS_8.rds")
# cmenas_2_WE <- readRDS("time_series_outputs/new/cmenas_I_9.rds")
# chondrus_2_WE <- readRDS("time_series_outputs/new/chondrus_I_8.rds")
# ascophyllum_2_WE <- readRDS("time_series_outputs/new/ascophyllum_GS_8.rds")
# 
# 
# # 2nd difference penalisation
# 
# mytilus_2 <- readRDS("time_series_outputs/output/mytilus_I_8.rds")
# semibalanus_2 <- readRDS("time_series_outputs/new/semibalanus_GS_8.rds")
# nucella_2 <- readRDS("time_series_outputs/new/nucella_I_8.rds")
# littorina_2 <- readRDS("time_series_outputs/new/littorina_I_10.rds")
# mastocarpus_2 <- readRDS("time_series_outputs/new/mastocarpus_GI_8.rds")
# fucus_2 <- readRDS("time_series_outputs/new/fucus_GS_8.rds")
# cmenas_2 <- readRDS("time_series_outputs/new/cmenas_I_9.rds")
# chondrus_2 <- readRDS("time_series_outputs/new/chondrus_I_8.rds")
# ascophyllum_2 <- readRDS("time_series_outputs/new/ascophyllum_GS_8.rds")
# 
# 
# # 3rd difference penalisation Compute Canada
# 
# mytilus_3_CC <- readRDS("time_series_outputs/new/mytilus_I_8.rds")
# semibalanus_3_CC <- readRDS("time_series_outputs/new/semibalanus_GS_8.rds")
# nucella_3_CC <- readRDS("time_series_outputs/new/nucella_I_8.rds")
# littorina_3_CC <- readRDS("time_series_outputs/new/littorina_I_8.rds")
# mastocarpus_3_CC <- readRDS("time_series_outputs/new/mastocarpus_GI_8.rds")
# fucus_3_CC <- readRDS("time_series_outputs/new/fucus_GS_8.rds")
# cmenas_3_CC <- readRDS("time_series_outputs/new/cmenas_I_9.rds")
# chondrus_3_CC <- readRDS("time_series_outputs/new/chondrus_I_8.rds")
# ascophyllum_3_CC <- readRDS("time_series_outputs/new/ascophyllum_GS_8.rds")
# 
# # 3rd difference local
# 
# mytilus_3_local <- readRDS("time_series_outputs/new/mytilus_I_9.rds")
# semibalanus_2_local <- readRDS("time_series_outputs/output/semibalanus_GS_6.rds")
# 
# 
# 
# # 3rd diff cc //
# 
# mytilus_3_paral <- readRDS("time_series_outputs/output/Fucked_up_server_parallel_run/mytilus_I_8.rds")


### 2nd difference no Wave exposure all I type

mytilus <- readRDS("time_series_outputs/new/mytilus_I_10_2.rds")
semibalanus <- readRDS("time_series_outputs/new/semibalanus_I_10.rds")
nucella <- readRDS("time_series_outputs/new/nucella_I_10.rds")
littorina <- readRDS("time_series_outputs/output/littorina_I_8.rds")
mastocarpus <- readRDS("time_series_outputs/new/mastocarpus_I_10.rds")
fucus <- readRDS("time_series_outputs/new/fucus_I_10.rds")
cmaenas <- readRDS("time_series_outputs/new/cmaenas_I_10.rds")
chondrus <- readRDS("time_series_outputs/new/chondrus_I_10.rds")
ascophyllum <- readRDS("time_series_outputs/new/ascophyllum_I_10.rds")

littorina <- readRDS("time_series_outputs/new/littorina_I_10.rds")
littorina <- readRDS("time_series_outputs/new/littorina_I_10_poisson.rds")

I <- mytilus #good!
I <- semibalanus #good!
I <- nucella #good!
I <- littorina # good with poisson distribution!!!!
I <- mastocarpus # good!
I <- fucus #good!
I <- cmaenas # converges well, but the fit is not very good. Essayer poisson - aide pas
I <- chondrus #good!
I <- ascophyllum #good!!


### Test convergence # # # # #  # #

I
# mcmc_plot(GS, type = "trace")
# mcmc_plot(GI, type = "trace")
mcmc_plot(I, type = "trace")

#DOES THE POSTERIOR INFORMATION HISTOGRAME HAVE ENOUGH INFO?

# mcmc_plot(GS, type = "hist")
# mcmc_plot(GI, type = "hist")
mcmc_plot(I, type = "hist")

# DO THE CHAIN EXHIBIT A STRONG DEGREE OF AUTO-CORRELATION?

# modelposterior_GS <- as.mcmc(GS)
# modelposterior_GI <- as.mcmc(GI)
modelposterior_I <- as.mcmc(I)

## look for temporal auto-cor. Does not show any.
# mcmc_plot(GS, type = "acf", pars = "zi_syear_1")
# mcmc_plot(GS, type = "acf", pars = "syear_1")
# mcmc_plot(GS, type = "acf", pars = "s(year)")
# mcmc_plot(GS, type = "acf", coef = "year")
# 
# mcmc_plot(GI, type = "acf", pars = "zi_syear_1") 
# mcmc_plot(GI, type = "acf", pars = "syear_1")
# mcmc_plot(GI, type = "acf", pars = "s(year)")
# mcmc_plot(GI, type = "acf", coef = "year")
# 

mcmc_plot(I, type = "acf", pars = "zi_syear:siteBabbsCove_1")
#mcmc_plot(I, type = "acf", pars = "zi_syear:siteMalagaCut_1")
mcmc_plot(I, type = "acf", pars = "zi_syear:siteNEAppledore_1") 
mcmc_plot(I, type = "acf", pars = "zi_syear:siteNWAppledore_1")
mcmc_plot(I, type = "acf", pars = "zi_syear:siteSEAppledore_1") 
mcmc_plot(I, type = "acf", pars = "zi_syear:siteSWAppledore_1")
mcmc_plot(I, type = "acf", pars = "syear:siteBabbsCove_1")
#mcmc_plot(I, type = "acf", pars = "syear:siteMalagaCut_1")
mcmc_plot(I, type = "acf", pars = "syear:siteNEAppledore_1") 
mcmc_plot(I, type = "acf", pars = "syear:siteNWAppledore_1")
mcmc_plot(I, type = "acf", pars = "syear:siteSEAppledore_1") 
mcmc_plot(I, type = "acf", pars = "syear:siteSWAppledore_1")


# mcmc_plot(GS, type = "acf", pars = "syear_1")
# mcmc_plot(GS, type = "acf", pars = "zi_syear_1")
# mcmc_plot(GS, type = "acf", pars = "t2yearsite_1")
# mcmc_plot(GS, type = "acf", pars = "t2yearsite_2")
# mcmc_plot(GS, type = "acf", pars = "t2yearsite_3")
# mcmc_plot(GS, type = "acf", pars = "zi_t2yearsite_1")
# mcmc_plot(GS, type = "acf", pars = "zi_t2yearsite_2")
# mcmc_plot(GS, type = "acf", pars = "zi_t2yearsite_3")
# mcmc_plot(GS, type = "acf", pars = "slevel_1")
# 
# 
# mcmc_plot(GI, type = "acf", pars = "syear_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syear_1")
# mcmc_plot(GI, type = "acf", pars = "syearsiteBabbsCove_1")
# mcmc_plot(GI, type = "acf", pars = "syearsiteMalagaCut_1")
# mcmc_plot(GI, type = "acf", pars = "syearsiteNEAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "syearsiteNWAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "syearsiteSEAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "syearsiteSWAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syearsiteBabbsCove_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syearsiteMalagaCut_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syearsiteNEAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syearsiteNWAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syearsiteSEAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "zi_syearsiteSWAppledore_1")
# mcmc_plot(GI, type = "acf", pars = "ssite_1")
# mcmc_plot(GI, type = "acf", pars = "zi_ssite_1")


## Plot posterior distribution
# 
# array_posterior_GS <- as.array(modelposterior_GS)
# array_posterior_GI <- as.array(modelposterior_GI)
array_posterior_I <- as.array(modelposterior_I)


# mcmc_areas(modelposterior_GS, pars = c("sds_syear_1", "sds_zi_syear_1"))
# mcmc_areas(modelposterior_GI, pars = c("sds_syear_1", "sds_zi_syear_1"))
mcmc_areas(modelposterior_I, pars = c("sds_syearsiteBabbsCove_1", "sds_zi_syearsiteBabbsCove_1",
                                      "sds_syearsiteNEAppledore_1", "sds_zi_syearsiteNEAppledore_1",
                                      "sds_syearsiteNWAppledore_1", "sds_zi_syearsiteNWAppledore_1",
                                      "sds_syearsiteSEAppledore_1", "sds_zi_syearsiteSEAppledore_1",
                                      "sds_syearsiteSWAppledore_1", "sds_zi_syearsiteSWAppledore_1"))

# mcmc_areas(modelposterior_GS, pars = c("sds_slevel_1"))
# mcmc_areas(modelposterior_GI, pars = c("sds_slevel_1"))
mcmc_areas(modelposterior_I, pars = c("sds_slevel_1"))

# mcmc_areas(modelposterior_GS, pars = c("sds_t2yearsite_1", "sds_zi_t2yearsite_1",
#                                        "sds_t2yearsite_2", "sds_zi_t2yearsite_2",
#                                        "sds_t2yearsite_3", "sds_zi_t2yearsite_3"))
# 
# mcmc_areas(modelposterior_GI, pars = c("sds_syearsiteBabbsCove_1", "sds_zi_syearsiteBabbsCove_1",
#                                        "sds_syearsiteMalagaCut_1", "sds_zi_syearsiteMalagaCut_1",
#                                        "sds_syearsiteNEAppledore_1", "sds_zi_syearsiteNEAppledore_1",
#                                        "sds_syearsiteNWAppledore_1", "sds_zi_syearsiteNWAppledore_1",
#                                        "sds_syearsiteSEAppledore_1", "sds_zi_syearsiteSEAppledore_1",
#                                        "sds_syearsiteSWAppledore_1", "sds_zi_syearsiteSWAppledore_1",
#                                        "sds_ssite_1", "sds_zi_ssite_1"))



## Check the fit
# 
# pp_check(GS, type = "ecdf_overlay")
# pp_check(GI, type = "ecdf_overlay")
pp_check(I, type = "ecdf_overlay")

# pp_check(GS, type = "ribbon_grouped", group = "site", nsamples = 25) 
# pp_check(GI, type = "ribbon_grouped", group = "site", nsamples = 25)
pp_check(I, type = "ribbon_grouped", group = "site", nsamples = 25) 


#pp_check(GS, type = "ribbon_grouped", group = "exposed", nsamples = 25) 
#pp_check(GI, type = "ribbon_grouped", group = "exposed", nsamples = 25)
#pp_check(I, type = "ribbon_grouped", group = "exposed", nsamples = 25) 


# Compare models with BIC

# GS$criteria$loo
# GI$criteria$loo
I$criteria$loo 

# GS$criteria$waic
# GI$criteria$waic
I$criteria$waic

mean(I$criteria$bayes_R2)

conditional_smooths(I)

