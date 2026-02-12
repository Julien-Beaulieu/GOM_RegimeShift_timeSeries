# This script is to fit multivariate time-series (GAM) on intertidal data and also 
# subtidal data from NOAA and compare the rates of changes.
# Contributor: Julien Beaulieu

# Packages

library(brms)

# load data

score_noaa <- read.csv("time_series_outputs/ordination_output/score_noaa_pcoa_site.csv")
cover_merge <- read.csv("time_series_outputs/ordination_output/score_PCoAsite_cover.csv")
count_merge <- read.csv("time_series_outputs/ordination_output/score_PCoAsite_count.csv")



cover_merge$INTERTIDAL_TRANSECT <- as.factor(cover_merge$INTERTIDAL_TRANSECT)
count_merge$INTERTIDAL_TRANSECT <- as.factor(count_merge$INTERTIDAL_TRANSECT)
score_noaa$year <- as.numeric(score_noaa$year)


# noaa model ==================================================================

#noaa_pcoa_mod <- bf( mvbind(MDS1,MDS2,MDS3) ~ s(year), family = gaussian())
#get_prior(noaa_pcoa_mod, data = score.pcoa.sites)

#prior_noaa_pcoa_mod <- c(prior(normal(0,1), class = "b", resp = ""),
#                         prior(normal(0,1), class = "b", resp ="MDS1"),
#                         prior(normal(0,1), class = "b", resp = "MDS2"),
#                         prior(normal(0,1), class = "b", resp = "MDS3"),
#                         prior(normal(0,1), class = "Intercept", resp = ""),
#                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS1"),
#                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS2"),
#                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS3"),
#                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS1"),
#                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS2"),
#                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS3"),
#                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS1"),
#                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS2"),
#                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS3"))

## fit !!!!


#fit <- brm(noaa_pcoa_mod,
#           prior = prior_noaa_pcoa_mod,
#           data = score_noaa,
#           warmup = 3000, iter = 7000,
#           thin = 20,
#           cores = 3, chains = 3,
#           control = list(adapt_delta = 0.99))

#saveRDS(fit, file = "time_series_outputs/noaa_pcoa_GS_3.rds")


# intertidal count =================================================================


GS <- bf( mvbind(MDS1,MDS2,MDS3) ~ 1 + s(YEAR,bs = "tp")+ s(LEVEL, bs = "tp") + t2(YEAR, SITE, bs = c("tp","re"), full = T) +(1|r|INTERTIDAL_TRANSECT),
          family = gaussian())
get_prior(GS, data = count_merge)


prior_pcoa_count <- c(prior(normal(0,1), class = "b", resp = ""),
                      prior(normal(0,1), class = "b", resp ="MDS1"),
                      prior(normal(0,1), class = "b", resp = "MDS2"),
                      prior(normal(0,1), class = "b", resp = "MDS3"),
                      prior(normal(0,1), class = "Intercept", resp = ""),
                      prior(student_t(3,-0.6,2.5), class = "Intercept", resp = "MDS1"),
                      prior(student_t(3,0.2,2.5), class = "Intercept", resp = "MDS2"),
                      prior(student_t(3,0.3,2.5), class = "Intercept", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS3"))
### FIT

fit <- brm(GS,
           prior = prior_pcoa_count,
           data = count_merge,
           warmup = 2000, iter = 7000,
           cores = 3, chains = 3,
           thin = 20,
           control = list(adapt_delta = 0.99))
saveRDS(fit, file = "time_series_outputs/output/count_pcoa_GS_10.rds")

##GI model (Group level trends different smoothness)

mod_GI <- bf(mvbind(MDS1,MDS2,MDS3) ~ 1+ s(YEAR,bs = "tp")+ s(YEAR, by = SITE, m=1, bs="tp") + s(LEVEL, bs = "tp") + s(SITE, bs = "re") +(1|r|INTERTIDAL_TRANSECT),
             family = gaussian())
get_prior(mod_GI, data = count_merge) # same as GS

fit <- brm(mod_GI,
           prior = prior_pcoa_count,
           data = count_merge,
           warmup = 2000, iter = 7000,
           cores = 3, chains = 3,
           thin = 20,
           control = list(adapt_delta = 0.99))
saveRDS(fit, file = "time_series_outputs/output/count_pcoa_GI_10.rds")

# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(mvbind(MDS1,MDS2,MDS3) ~ 1 + s(YEAR, by = SITE, bs="tp") + s(LEVEL, bs = "tp") + (1|r|INTERTIDAL_TRANSECT),
            family = gaussian())
get_prior(mod_I, data = count_merge) # same as GS


fit <- brm(mod_I,
           prior = prior_pcoa_count,
           data = count_merge,
           warmup = 2000, iter = 7000,
           cores = 3, chains = 3,
           thin = 20,
           control = list(adapt_delta = 0.99))

saveRDS(fit, file = "time_series_outputs/output/count_pcoa_I_10.rds")

## Model with best mix of types =======================

MDS1 <- bf(MDS1 ~ 1 + s(YEAR, bs = "tp",m=3) + s(LEVEL, bs = "tp") + 
             t2(YEAR, SITE, bs = c("tp", "re"), full = T) + 
             (1 | r | INTERTIDAL_TRANSECT),family = gaussian())
MDS2 <- bf(MDS2 ~ 1 + s(YEAR, by = SITE, bs = "tp",m=3) + s(LEVEL, bs = "tp") +
             (1 | r | INTERTIDAL_TRANSECT),family = gaussian())
MDS3 <- bf(MDS3 ~ 1 + s(YEAR, bs = "tp",m=3) + s(LEVEL, bs = "tp") + 
             t2(YEAR, SITE, bs = c("tp", "re"), full = T) + 
             (1 | r | INTERTIDAL_TRANSECT),family = gaussian())

mod_mix <- mvbrmsformula(MDS1,MDS2,MDS3)
get_prior(mod_mix, data = count_merge)

prior_pcoa_count <- c(prior(normal(0,1), class = "b", resp = ""),
                      prior(normal(0,1), class = "b", resp ="MDS1"),
                      prior(normal(0,1), class = "b", resp = "MDS2"),
                      prior(normal(0,1), class = "b", resp = "MDS3"),
                      prior(normal(0,1), class = "Intercept", resp = ""),
                      prior(student_t(3,-0.3,2.5), class = "Intercept", resp = "MDS1"),
                      prior(student_t(3,0.2,2.5), class = "Intercept", resp = "MDS2"),
                      prior(student_t(3,0.3,2.5), class = "Intercept", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sd", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sd", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sd", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS3"))
## FIt

fit <- brm(mod_mix,
           prior = prior_pcoa_count,
           data = count_merge,
           warmup = 2000, iter = 7000,
           cores = 4, chains = 4,
           thin = 20,
           control = list(adapt_delta = 0.99))
fit <- add_criterion(fit, c("loo","bayes_R2","loo_R2"))

saveRDS(fit, file = "time_series_outputs/new/count_pcoa_mix_8.rds")



# Intertidal % cover ===========================================================

# GS <- bf( mvbind(MDS1,MDS2,MDS3) ~ 1 + s(YEAR,bs = "tp",m=3)+ s(LEVEL, bs = "tp") + t2(YEAR, SITE, bs = c("tp","re"), full = T) +(1|r|INTERTIDAL_TRANSECT),
#           family = gaussian())
#get_prior(GS, data = cover_merge)


prior_pcoa_cover <- c(prior(normal(0,1), class = "b", resp = ""),
                      prior(normal(0,1), class = "b", resp ="MDS1"),
                      prior(normal(0,1), class = "b", resp = "MDS2"),
                      prior(normal(0,1), class = "b", resp = "MDS3"),
                      prior(normal(0,1), class = "Intercept", resp = ""),
                      prior(student_t(3,0.4,2.5), class = "Intercept", resp = "MDS1"),
                      prior(student_t(3,0.1,2.5), class = "Intercept", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sds", resp = "MDS3"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS1"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS2"),
                      prior(student_t(3,0,2.5), class = "sigma", resp = "MDS3"))

# FIT

# fit <- brm(GS,
#            prior = prior_pcoa_cover,
#            data = cover_merge,
#            warmup = 2000, iter = 7000,
#            cores = 4, chains = 4,
#            thin = 20,
#            control = list(adapt_delta = 0.99))
# fit <- add_criterion(fit, c("loo","bayes_R2","loo_R2"))
# saveRDS(fit, file = "time_series_outputs/output/cover_pcoa_GS_8.rds")


##GI model (Group level trends different smoothness)

mod_GI <- bf(mvbind(MDS1,MDS2,MDS3) ~ 1 + s(YEAR,bs = "tp")+ s(YEAR, by = SITE, m=1, bs="tp") + s(LEVEL, bs = "tp") + s(SITE, bs = "re") +(1|SITE/INTERTIDAL_TRANSECT),
             family = gaussian())
get_prior(mod_GI, data = cover_merge) # same as GS

fit <- brm(mod_GI,
           prior = prior_pcoa_cover,
           data = cover_merge,
           warmup = 2000, iter = 7000,
           cores = 3, chains = 3,
           thin = 20,
           control = list(adapt_delta = 0.99))

saveRDS(fit, file = "time_series_outputs/output/cover_pcoa_GI_10.rds")


# # I model (no shared trends and sifferent smoothness)

mod_I <- bf(mvbind(MDS1,MDS2,MDS3) ~ 1 + exposed + s(YEAR, by = SITE, bs = "tp", m = 2) + s(LEVEL, bs = "tp") + (1|SITE/INTERTIDAL_TRANSECT),
            family = gaussian())
#get_prior(mod_I, data = cover_merge) # same as GS


fit <- brm(mod_I,
           prior = prior_pcoa_cover,
           data = cover_merge,
           warmup = 2000, iter = 7000,
           cores = 3, chains = 3,
           thin = 20,
           control = list(adapt_delta = 0.99))

saveRDS(fit, file = "time_series_outputs/output/cover_pcoa_I_10.rds")

