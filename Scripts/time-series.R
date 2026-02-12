# Time series
# This script is to fit time series on the principal organism from the SEMs.
# Contributor: Julien Beaulieu

# 1- Load data
# 2- Load packages 
# 3- Mussels
# 4- Barnacle
# 5- fucus
# 6- ascophylum
# 7- chondrus
# 8- littorina
# 9- mastocarpus
# 10- green crab
# 11- nucella


# 1: load packages #################################################


library(brms)
library(schoenberg)
library(readr)
library(coda)
library(boot)
library(dplyr)
library(tidyverse)

# set cores for running on server
#n_cores <- parallel::detectCores()
n_chains <- 4
#n_threads <- n_cores/n_chains
n_cores = 4

# 2: Load data 

dat <- read_csv("data/full_dataset_for_SEMs_groupcentered.csv", guess_max = 3000) #data produced by SEM_dataset_assembly.R
dat$site <- as.character(dat$site)
dat$intertidal_transect <- as.factor(dat$intertidal_transect)
dat <- filter(dat, site != "Malaga Cut") #remove malaga cut
dat <- dat[,-1]


#load raw intertidal data to get isopoda and amphipoda
intertidal_full <- read_rds("data/combined_intertidal_abundance.RDS")

amph_iso_dat <- intertidal_full %>% filter(ORGANISM == "Amphipoda"|ORGANISM == "Isopods") %>% 
  select(c(SITE,INTERTIDAL_TRANSECT,YEAR,LEVEL,REPLICATE,ORGANISM,VALUE)) %>%
  unite(ID, c(SITE,INTERTIDAL_TRANSECT,YEAR,LEVEL,ORGANISM,REPLICATE), sep = "/",remove = F) %>%
  distinct(ID, .keep_all = TRUE) %>% 
  separate(ID, c("site","intertidal_transect","year","level","specie","replicate"), sep = "/") %>% 
  dplyr::select(site,intertidal_transect,year,level,specie,replicate,VALUE) %>% 
  pivot_wider(names_from = specie, values_from = VALUE, values_fill = 0) %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(level = as.numeric(level)) %>% 
  filter(site != "Malaga Cut")

# Add rows for sum littorina and fucus 

#dat$sum_littorina <- rowSums(dat[,c("littorina_littorea", "littorina_obtusata", "littorina_saxatilis")])
#dat$sum_fucus <- rowSums(dat[,c("fucus_vesiculosus","fucus_distichus")])
#dat$sum_fucus <- as.integer(dat$sum_fucus)
dat$littorina <- as.integer(dat$littorina)


### count the numger of transect each year #######
y <- unique(dat$year)
for (i in c(1:length(y))) {
  dat2 <- dat %>% filter(year == y[i])
  nb_transect[i] <- length(unique(dat2$intertidal_transect))
  print(length(unique(dat2$intertidal_transect)))
}

# 3: mytilus_edulis######################################################
##________________________________________________________________________

#GS model (group level trends and similar smoothness)

#mod_GS <- bf(mytilus_edulis ~ 1 + exposed + s(year,bs = "tp")+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year, bs = "tp") + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
#             family = zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"))
#get_prior(mod_GS, data = dat) #same as BM

prior_mytilus <- c(prior(student_t(3, 0, 2.5), class = "Intercept"),
                 prior(normal(0, 1), class = "b"),
                 prior(gamma(0.01, 0.01), class = "phi"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                 prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                 prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                 prior(student_t(3,0,2.5), class = "sd", dpar = "zi"), # SD of the random effect -- use exponential one
                 prior(student_t(3,0,2.5), class = "sds", dpar = "zi"),
                 prior(logistic(0,1), class = "Intercept", dpar = "zi"),
                 prior(normal(0, 1), class = "b", dpar = "zi")
)


# FIT

#fit <- brm(mod_GS,
#           prior = prior_mytilus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/mytilus_GS_3.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(mytilus_edulis ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"))
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_mytilus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/mytilus_GI_3.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(mytilus_edulis ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"))
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_mytilus,
           data = dat,
           warmup = 2000, iter = 5000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/mytilus_I_10_2.rds")


print("mytilus done")



# 4: Semibalanus_balanoides######################################################
##________________________________________________________________________


#GS model (group level trends and similar smoothness)

mod_GS <- bf(semibalanus_balanoides ~ 1 + s(year,bs = "tp", m = 2)+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
             zi ~ 1 + s(year, bs = "tp", m = 2) + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
             family = zero_inflated_beta())
#get_prior(mod_GS, data = dat) #same as BM

prior_semibalanus <- c(prior(student_t(3, 0, 2.5), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(gamma(0.01, 0.01), class = "phi"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                   prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                   prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                   prior(student_t(3,0,2.5), class = "sd", dpar = "zi"), # SD of the random effect -- use exponential one
                   prior(student_t(3,0,2.5), class = "sds", dpar = "zi"),
                   prior(logistic(0,1), class = "Intercept", dpar = "zi"),
                   prior(normal(0, 1), class = "b", dpar = "zi")
)


# FIT

fit <- brm(mod_GS,
           prior = prior_semibalanus,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/semibalanus_GS_10.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(semibalanus_balanoides ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = zero_inflated_beta())
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_semibalanus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/semibalanus_GI_5.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(semibalanus_balanoides ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = zero_inflated_beta())
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_semibalanus,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
saveRDS(fit, file = "time_series_outputs/new/semibalanus_I_10.rds")


print("semibalanus done")


# 5: sum_fucus ##########################################################
#______________________________________________________________________

#GS model (group level trends and similar smoothness)

mod_GS <- bf(fucus ~ 1 + s(year,bs = "tp",m = 2)+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
             hu ~ 1 + s(year, bs = "tp",m = 2) + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
             family = hurdle_gamma(link = "log"))
get_prior(mod_GS, data = dat) #same as BM

prior_fucus <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                       prior(normal(0, 1), class = "b"),
                       prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                       prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                       prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                       prior(student_t(3,0,2.5), class = "sd", dpar = "hu"), # SD of the random effect -- use exponential one
                       prior(student_t(3,0,2.5), class = "sds", dpar = "hu"),
                       prior(logistic(0,1), class = "Intercept", dpar = "hu"),
                       prior(normal(0, 1), class = "b", dpar = "hu")
)


# FIT

fit <- brm(mod_GS,
           prior = prior_fucus,
           data = dat,
           warmup = 2000, iter = 5000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/fucus_GS_10.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(sum_fucus ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             hu ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = hurdle_gamma(link = "log"))
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_fucus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/fucus_GI_3.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(fucus ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            hu ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = hurdle_gamma(link = "log"))
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_fucus,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
saveRDS(fit, file = "time_series_outputs/new/fucus_I_10.rds")


print("fucus done")


# 6: ascophylum_nodusum #####################################################
#_____________________________________________________________________________

#GS model (group level trends and similar smoothness)

mod_GS <- bf(ascophyllum_nodosum ~ 1 + s(year,bs = "tp",m = 2) + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
             hu ~ 1 + s(year, bs = "tp",m = 2) + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
             family = hurdle_gamma(link = "log"))
#get_prior(mod_GS, data = dat) #same as BM

prior_ascophyllum <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                       prior(normal(0, 1), class = "b"),
                       prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                       prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                       prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                       prior(student_t(3,0,2.5), class = "sd", dpar = "hu"), # SD of the random effect -- use exponential one
                       prior(student_t(3,0,2.5), class = "sds", dpar = "hu"),
                       prior(logistic(0,1), class = "Intercept", dpar = "hu"),
                       prior(normal(0, 1), class = "b", dpar = "hu")
)


# FIT

fit <- brm(mod_GS,
           prior = prior_ascophyllum,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/ascophyllum_GS_10.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(ascophyllum_nodosum ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             hu ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = hurdle_gamma(link = "log"))
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_ascophyllum,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/ascophyllum_GI_3.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(ascophyllum_nodosum ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            hu ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = hurdle_gamma(link = "log"))
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_ascophyllum,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
saveRDS(fit, file = "time_series_outputs/new/ascophyllum_I_10.rds")


print("ascophyllum done")

# 7: chondrus ############################################################
#________________________________________________________________________

#GS model (group level trends and similar smoothness)

#mod_GS <- bf(chondrus_crispus ~ 1 + exposed + s(year,bs = "tp")+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
#             hu ~ 1 + exposed + s(year, bs = "tp") + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
#             family = hurdle_gamma(link = "log"))
#get_prior(mod_GS, data = dat) #same as BM

prior_chondrus <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                       prior(normal(0, 1), class = "b"),
                       prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                       prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                       prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                       prior(student_t(3,0,2.5), class = "sd", dpar = "hu"), # SD of the random effect -- use exponential one
                       prior(student_t(3,0,2.5), class = "sds", dpar = "hu"),
                       prior(logistic(0,1), class = "Intercept", dpar = "hu"),
                       prior(normal(0, 1), class = "b", dpar = "hu")
)


# FIT

#fit <- brm(mod_GS,
#           prior = prior_chondrus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/chondrus_GS_3.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(chondrus_crispus ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             hu ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = hurdle_gamma(link = "log"))
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_chondrus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2"))
#saveRDS(fit, file = "output/chondrus_GI_3.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(chondrus_crispus ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            hu ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = hurdle_gamma(link = "log"))
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_chondrus,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/chondrus_I_10.rds")

print("chondrus done")

# 8: littorina ##################################################################
#____________________________________________________________________________

#GS model (group level trends and similar smoothness)

#mod_GS <- bf(littorina ~ 1 + exposed + s(year,bs = "tp")+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year, bs = "tp") + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
#             family = zero_inflated_negbinomial())
#get_prior(mod_GS, data = dat) #same as BM

prior_littorina <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                    prior(normal(0,1), class = "b"),
                    #prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                    prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                    prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                    prior(student_t(3,0,2.5), class = "sd", dpar = "zi"), # SD of the random effect -- use exponential one
                    prior(student_t(3,0,2.5), class = "sds", dpar = "zi"),
                    prior(logistic(0,1), class = "Intercept", dpar = "zi"),
                    prior(normal(0, 1), class = "b", dpar = "zi")
)


### FIT

#fit <- brm(mod_GS,
#           prior = prior_littorina,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/littorina_GS_5.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(littorina ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = zero_inflated_negbinomial())
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_littorina,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/littorina_GI_5.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(littorina ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = zero_inflated_poisson())
get_prior(mod_I, data = dat2) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_littorina,
           data = dat,
           warmup = 2000, iter = 9000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/littorina_I_10_poisson.rds")


print("littorina done")


# 9: mastocarpus##############################################################
#___________________________________________________________________________

#GS model (group level trends and similar smoothness)

#mod_GS <- bf(mastocarpus_stellatus ~ 1 + exposed + s(year,bs = "tp")+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
#             hu ~ 1 + exposed + s(year, bs = "tp") + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
#             family = hurdle_gamma(link = "log"))
#get_prior(mod_GS, data = dat) #same as BM

prior_mastocarpus <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                    prior(normal(0, 1), class = "b"),
                    prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                    prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                    prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                    prior(student_t(3,0,2.5), class = "sd", dpar = "hu"), # SD of the random effect -- use exponential one
                    prior(student_t(3,0,2.5), class = "sds", dpar = "hu"),
                    prior(logistic(0,1), class = "Intercept", dpar = "hu"),
                    prior(normal(0, 1), class = "b", dpar = "hu")
)


# FIT

#fit <- brm(mod_GS,
#           prior = prior_mastocarpus,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/mastocarpus_GS_3.rds")


#GI model (Group level trends different smoothness)

mod_GI <- bf(mastocarpus_stellatus ~ 1+ s(year,bs = "tp",m = 2)+ s(year, by = site, bs="tp",m = 2) + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
             hu ~ 1 + s(year,bs = "tp",m = 2)+ s(year, by = site, bs="tp",m = 2) + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
             family = hurdle_gamma(link = "log"))
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
fit <- brm(mod_GI,
           prior = prior_mastocarpus,
           data = dat,
           warmup = 2000, iter = 6000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.999)) 
fit <- add_criterion(fit, c("waic","loo","loo_R2","bayes_R2"))
saveRDS(fit, file = "time_series_outputs/new/mastocarpus_GI_10.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(mastocarpus_stellatus ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            hu ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = hurdle_gamma(link = "log"))
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_mastocarpus,
           data = dat,
           warmup = 2000, iter = 7000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
saveRDS(fit, file = "time_series_outputs/new/mastocarpus_I_10.rds")


print("mastocarpus done")

# 10: green crab ###########################################################
#___________________________________________________________________________

#GS model (group level trends and similar smoothness)

#mod_GS <- bf(carcinus_maenas ~ 1 + exposed + s(year,bs = "tp")+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year, bs = "tp") + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
#             family = zero_inflated_negbinomial())
#get_prior(mod_GS, data = dat) #same as BM

prior_cmenas <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     #prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                     prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                     prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                     prior(student_t(3,0,2.5), class = "sd", dpar = "zi"), # SD of the random effect -- use exponential one
                     prior(student_t(3,0,2.5), class = "sds", dpar = "zi"),
                     prior(logistic(0,1), class = "Intercept", dpar = "zi"),
                     prior(normal(0, 1), class = "b", dpar = "zi")
)


# FIT

#fit <- brm(mod_GS,
#           prior = prior_cmenas,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/cmenas_GS_3.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(carcinus_maenas ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=1, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = zero_inflated_negbinomial())
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_cmenas,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/cmenas_GI_3.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(carcinus_maenas ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp", m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            #family = zero_inflated_negbinomial())
            family = zero_inflated_poisson())
get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_cmenas,
           data = dat,
           warmup = 3000, iter = 9000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/cmaenas_I_10_poisson.rds")
fit <- readRDS("time_series_outputs/new/cmaenas_I_10_poisson.rds")

print("cmenas done")
 

# 11: Nucella ##############################################################
#___________________________________________________________________________

#GS model (group level trends and similar smoothness)

#mod_GS <- bf(nucella_lapillus ~ 1 + exposed + s(year,bs = "tp")+ s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year, bs = "tp") + s(level, bs = "tp") + t2(year, site, bs = c("tp","re"), full = T) + (1|r|intertidal_transect),
#             family = zero_inflated_negbinomial())
#get_prior(mod_GS, data = dat) #same as BM

prior_nucella <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                  prior(normal(0, 1), class = "b"),
                  prior(gamma(0.01, 0.01), class = "shape"), # no need to duplicate for zi because only one estimate of phi = precision parameter: if very certain about the mean phi is high and if low phi is low
                  prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                  prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                  prior(student_t(3,0,2.5), class = "sd", dpar = "zi"), # SD of the random effect -- use exponential one
                  prior(student_t(3,0,2.5), class = "sds", dpar = "zi"),
                  prior(logistic(0,1), class = "Intercept", dpar = "zi"),
                  prior(normal(0, 1), class = "b", dpar = "zi")
)


# FIT

#fit <- brm(mod_GS,
#           prior = prior_nucella,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/nucella_GS_5.rds")


#GI model (Group level trends different smoothness)

#mod_GI <- bf(nucella_lapillus ~ 1+ exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             zi ~ 1 + exposed + s(year,bs = "tp")+ s(year, by = site, m=2, bs="tp") + s(level, bs = "tp") + s(site, bs = "re") +(1|r|intertidal_transect),
#             family = zero_inflated_negbinomial())
#get_prior(mod_GI, data = dat) # same as GS


# FIT # # # # # # # # # # # # # # # 
#fit <- brm(mod_GI,
#           prior = prior_nucella,
#           data = dat,
#           warmup = 2000, iter = 7000,
#           cores = n_cores, chains = n_chains,
#           backend = "cmdstanr", threads = threading(n_threads),
#           seed = 17,
#           control = list(adapt_delta = 0.99)) 
#fit <- add_criterion(fit, c("loo","waic","loo_R2","bayes_R2"))
#saveRDS(fit, file = "output/nucella_GI_5.rds")



# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(nucella_lapillus ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = zero_inflated_negbinomial())
#get_prior(mod_I, data = dat) # same as GS


# FIT
fit <- brm(mod_I,
           prior = prior_nucella,
           data = dat,
           warmup = 3000, iter = 9000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/nucella_I_10.rds")


print("nucella done")


#### Amphipoda #################################################
####_____________________________________________________________

mod_I <- bf(Amphipoda ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = zero_inflated_poisson())
get_prior(mod_I, data = amph_iso_dat)

prior_amph <- c(prior(student_t(3, -2.3, 2.5), class = "Intercept"),
                prior(normal(0, 1), class = "b"),
                prior(student_t(3,0,2.5), class = "sd"), # SD of the random effect -- use exponential one
                prior(student_t(3,0,2.5), class = "sds"), # number of points in non linear interpolation.. number of exponents?
                prior(student_t(3,0,2.5), class = "sd", dpar = "zi"), # SD of the random effect -- use exponential one
                prior(student_t(3,0,2.5), class = "sds", dpar = "zi"),
                prior(logistic(0,1), class = "Intercept", dpar = "zi"),
                prior(normal(0, 1), class = "b", dpar = "zi")
)

fit <- brm(mod_I,
           data = amph_iso_dat,
           prior = prior_amph,
           warmup = 3000, iter = 9000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.9)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/Amphipoda_I_3.rds")

print("amphipoda done")

#### Amphipoda #################################################
####_____________________________________________________________

mod_I <- bf(Isopods ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            zi ~ 1 + s(year, by = site, bs = "tp",m = 2) + s(level, bs = "tp") + (1|r|intertidal_transect),
            family = zero_inflated_poisson())
get_prior(mod_I, data = amph_iso_dat)


fit <- brm(mod_I,
           data = amph_iso_dat,
           prior = prior_amph,
           warmup = 3000, iter = 9000,
           cores = n_cores, chains = n_chains,
           #backend = "cmdstanr", threads = threading(n_threads),
           seed = 17,
           control = list(adapt_delta = 0.99)) 
fit <- add_criterion(fit, c("loo_R2","bayes_R2","loo","waic"))
saveRDS(fit, file = "time_series_outputs/new/Isopods_I_2.rds")


print("finish")
