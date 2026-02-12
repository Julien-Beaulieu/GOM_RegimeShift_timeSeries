# This script is to calculate the FIRST derivative using forward difference.
# This script is to calculate the FIRST derivative to know if there is a change
# in intertidal communities when there is a regime shift. The script also produces
# a table with the 1st deriv values during subtidal regime shifts and produce 1st
# deriv figures.
# Contributor: Julien Beaulieu


#library
library(brms)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(bayestestR)

# load models results
# cover_GS <- readRDS("time_series_outputs/new/cover_pcoa_GS_9.rds")
# count_4 <- readRDS("time_series_outputs/new/count_pcoa_mix_9.rds")
count <- readRDS("time_series_outputs/output/count_pcoa_GI_10.rds")
cover <- readRDS("time_series_outputs/output/cover_pcoa_GI_10.rds")
noaa <- readRDS("time_series_outputs/output/noaa_pcoa_GS_3.rds")


mytilus <- readRDS("time_series_outputs/new/final_UV_TS/mytilus_I_10_2.rds")
semibalanus <- readRDS("time_series_outputs/new/final_UV_TS/semibalanus_I_10.rds")
nucella <- readRDS("time_series_outputs/new/final_UV_TS/nucella_I_10.rds") 
littorina <- readRDS("time_series_outputs/new/final_UV_TS/littorina_I_10_poisson.rds") 
mastocarpus <- readRDS("time_series_outputs/new/final_UV_TS/mastocarpus_I_10.rds")  
fucus <- readRDS("time_series_outputs/new/final_UV_TS/fucus_I_10.rds")
cmenas <- readRDS("time_series_outputs/new/final_UV_TS/cmaenas_I_10.rds") # convergence problem
chondrus <- readRDS("time_series_outputs/new/final_UV_TS/chondrus_I_10.rds")
ascophyllum <- readRDS("time_series_outputs/new/final_UV_TS/ascophyllum_I_10.rds")
amphipoda <- readRDS("time_series_outputs/new/final_UV_TS/Amphipoda_I_3.rds")
isopoda <- readRDS("time_series_outputs/new/final_UV_TS/Isopods_I_2.rds")


# calculate marginal pp for multivariate models ###########################
msms_cover <- conditional_smooths(cover, spaghetti = T)
msms_count <- conditional_smooths(count, spaghetti = T)


## Extract info from conditional_smooth========================================

### cover ##_______________________________________________________________________________
#extract data of the draws (spaghetti)

### cover global trend ###############################################


f <- attributes(msms_cover$`mu_MDS1: s(YEAR,bs="tp")`)$spaghetti
f <- f[,c("YEAR", "estimate__", "sample__")]
colnames(f)[1] <- "year"

global_cover_MDS1 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
global_cover_MDS1[is.nan(global_cover_MDS1)] <- 0


f <- attributes(msms_cover$`mu_MDS2: s(YEAR,bs="tp")`)$spaghetti
f <- f[,c("YEAR", "estimate__", "sample__")]
colnames(f)[1] <- "year"
global_cover_MDS2 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
global_cover_MDS2[is.nan(global_cover_MDS2)] <- 0

f <- attributes(msms_cover$`mu_MDS3: s(YEAR,bs="tp`)$spaghetti
f <- f[,c("YEAR", "estimate__", "sample__")]
colnames(f)[1] <- "year"
global_cover_MDS3 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
global_cover_MDS3[is.nan(global_cover_MDS3)] <- 0


### count global trend ############################################


f <- attributes(msms_count$`mu_MDS1: s(YEAR,bs="tp",m=3)`)$spaghetti
f <- f[,c("YEAR", "estimate__", "sample__")]
colnames(f)[1] <- "year"

global_count_MDS1 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
global_count_MDS1[is.nan(global_count_MDS1)] <- 0


f <- attributes(msms_count$`mu_MDS2: s(YEAR,bs="tp",m=3)`)$spaghetti
f <- f[,c("YEAR", "estimate__", "sample__")]
colnames(f)[1] <- "year"
global_count_MDS2 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
global_count_MDS2[is.nan(global_count_MDS2)] <- 0

f <- attributes(msms_count$`mu_MDS3: s(YEAR,bs="tp",m=3)`)$spaghetti
f <- f[,c("YEAR", "estimate__", "sample__")]
colnames(f)[1] <- "year"
global_count_MDS3 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
global_count_MDS3[is.nan(global_count_MDS3)] <- 0


### cover site specific #####################################

### Function to extract pp for univariate models

source("Scripts/pp_extract_function.R")

site_cover_MDS1 <- attributes(msms_cover$`mu_MDS1: s(YEAR,by=SITE,m=1,bs="tp")`)$spaghetti
pp_site_cover_MDS1 <- extract_site_MV(site_cover_MDS1)
pp_site_cover_MDS1_BC <- pp_site_cover_MDS1[[1]]
pp_site_cover_MDS1_NE <- pp_site_cover_MDS1[[2]]
pp_site_cover_MDS1_NW <- pp_site_cover_MDS1[[3]]
pp_site_cover_MDS1_SE <- pp_site_cover_MDS1[[4]]
pp_site_cover_MDS1_SW <- pp_site_cover_MDS1[[5]]

site_cover_MDS2 <- attributes(msms_cover$`mu_MDS2: s(YEAR,by=SITE,m=1,bs="tp")`)$spaghetti
pp_site_cover_MDS2 <- extract_site_MV(site_cover_MDS2)
pp_site_cover_MDS2_BC <- pp_site_cover_MDS2[[1]]
pp_site_cover_MDS2_NE <- pp_site_cover_MDS2[[2]]
pp_site_cover_MDS2_NW <- pp_site_cover_MDS2[[3]]
pp_site_cover_MDS2_SE <- pp_site_cover_MDS2[[4]]
pp_site_cover_MDS2_SW <- pp_site_cover_MDS2[[5]]

site_cover_MDS3 <- attributes(msms_cover$`mu_MDS3: s(YEAR,by=SITE,m=1,bs="tp")`)$spaghetti
pp_site_cover_MDS3 <- extract_site_MV(site_cover_MDS3)
pp_site_cover_MDS3_BC <- pp_site_cover_MDS3[[1]]
pp_site_cover_MDS3_NE <- pp_site_cover_MDS3[[2]]
pp_site_cover_MDS3_NW <- pp_site_cover_MDS3[[3]]
pp_site_cover_MDS3_SE <- pp_site_cover_MDS3[[4]]
pp_site_cover_MDS3_SW <- pp_site_cover_MDS3[[5]]


### Count site specific ####################################################

site_count_MDS1 <- attributes(msms_count$`mu_MDS1: s(YEAR,by=SITE,m=1,bs="tp")`)$spaghetti
pp_site_count_MDS1 <- extract_site_MV(site_count_MDS1)
pp_site_count_MDS1_BC <- pp_site_count_MDS1[[1]]
pp_site_count_MDS1_NE <- pp_site_count_MDS1[[2]]
pp_site_count_MDS1_NW <- pp_site_count_MDS1[[3]]
pp_site_count_MDS1_SE <- pp_site_count_MDS1[[4]]
pp_site_count_MDS1_SW <- pp_site_count_MDS1[[5]]


site_count_MDS2 <- attributes(msms_count$`mu_MDS2: s(YEAR,by=SITE,m=1,bs="tp")`)$spaghetti
pp_site_count_MDS2 <- extract_site_MV(site_count_MDS2)
pp_site_count_MDS2_BC <- pp_site_count_MDS2[[1]]
pp_site_count_MDS2_NE <- pp_site_count_MDS2[[2]]
pp_site_count_MDS2_NW <- pp_site_count_MDS2[[3]]
pp_site_count_MDS2_SE <- pp_site_count_MDS2[[4]]
pp_site_count_MDS2_SW <- pp_site_count_MDS2[[5]]

site_count_MDS3 <- attributes(msms_count$`mu_MDS3: s(YEAR,by=SITE,m=1,bs="tp")`)$spaghetti
pp_site_count_MDS3 <- extract_site_MV(site_count_MDS2)
pp_site_count_MDS3_BC <- pp_site_count_MDS3[[1]]
pp_site_count_MDS3_NE <- pp_site_count_MDS3[[2]]
pp_site_count_MDS3_NW <- pp_site_count_MDS3[[3]]
pp_site_count_MDS3_SE <- pp_site_count_MDS3[[4]]
pp_site_count_MDS3_SW <- pp_site_count_MDS3[[5]]


### Extract posterior prediction from all UV models#############################

# use the function for each specie 

pp_ascophylum <- extract(ascophyllum, model_type = "I_hu")
pp_ascophylum_BC <- pp_ascophylum[[1]]
pp_ascophylum_NE <- pp_ascophylum[[2]]
pp_ascophylum_NW <- pp_ascophylum[[3]]
pp_ascophylum_SE <- pp_ascophylum[[4]]
pp_ascophylum_SW <- pp_ascophylum[[5]]
pp_ascophylum_BC_zi <- pp_ascophylum[[6]]
pp_ascophylum_NE_zi <- pp_ascophylum[[7]]
pp_ascophylum_NW_zi <- pp_ascophylum[[8]]
pp_ascophylum_SE_zi <- pp_ascophylum[[9]]
pp_ascophylum_SW_zi <- pp_ascophylum[[10]]

pp_chondrus <- extract(chondrus, model_type = "I_hu")
pp_chondrus_BC <- pp_chondrus[[1]]
pp_chondrus_NE <- pp_chondrus[[2]]
pp_chondrus_NW <- pp_chondrus[[3]]
pp_chondrus_SE <- pp_chondrus[[4]]
pp_chondrus_SW <- pp_chondrus[[5]]
pp_chondrus_BC_zi <- pp_chondrus[[6]]
pp_chondrus_NE_zi <- pp_chondrus[[7]]
pp_chondrus_NW_zi <- pp_chondrus[[8]]
pp_chondrus_SE_zi <- pp_chondrus[[9]]
pp_chondrus_SW_zi <- pp_chondrus[[10]]

pp_cmaenas <- extract(cmenas, model_type = "I_zi")
pp_cmaenas_BC <- pp_cmaenas[[1]]
pp_cmaenas_NE <- pp_cmaenas[[2]]
pp_cmaenas_NW <- pp_cmaenas[[3]]
pp_cmaenas_SE <- pp_cmaenas[[4]]
pp_cmaenas_SW <- pp_cmaenas[[5]]
pp_cmaenas_BC_zi <- pp_cmaenas[[6]]
pp_cmaenas_NE_zi <- pp_cmaenas[[7]]
pp_cmaenas_NW_zi <- pp_cmaenas[[8]]
pp_cmaenas_SE_zi <- pp_cmaenas[[9]]
pp_cmaenas_SW_zi <- pp_cmaenas[[10]]

pp_fucus <- extract(fucus, model_type = "I_hu")
pp_fucus_BC <- pp_fucus[[1]]
pp_fucus_NE <- pp_fucus[[2]]
pp_fucus_NW <- pp_fucus[[3]]
pp_fucus_SE <- pp_fucus[[4]]
pp_fucus_SW <- pp_fucus[[5]]
pp_fucus_BC_zi <- pp_fucus[[6]]
pp_fucus_NE_zi <- pp_fucus[[7]]
pp_fucus_NW_zi <- pp_fucus[[8]]
pp_fucus_SE_zi <- pp_fucus[[9]]
pp_fucus_SW_zi <- pp_fucus[[10]]

pp_littorina <- extract(littorina, model_type = "I_zi")
pp_littorina_BC <- pp_littorina[[1]]
pp_littorina_NE <- pp_littorina[[2]]
pp_littorina_NW <- pp_littorina[[3]]
pp_littorina_SE <- pp_littorina[[4]]
pp_littorina_SW <- pp_littorina[[5]]
pp_littorina_BC_zi <- pp_littorina[[6]]
pp_littorina_NE_zi <- pp_littorina[[7]]
pp_littorina_NW_zi <- pp_littorina[[8]]
pp_littorina_SE_zi <- pp_littorina[[9]]
pp_littorina_SW_zi <- pp_littorina[[10]]

pp_mastocarpus <- extract(mastocarpus, model_type = "I_hu")
pp_mastocarpus_BC <- pp_mastocarpus[[1]]
pp_mastocarpus_NE <- pp_mastocarpus[[2]]
pp_mastocarpus_NW <- pp_mastocarpus[[3]]
pp_mastocarpus_SE <- pp_mastocarpus[[4]]
pp_mastocarpus_SW <- pp_mastocarpus[[5]]
pp_mastocarpus_BC_zi <- pp_mastocarpus[[6]]
pp_mastocarpus_NE_zi <- pp_mastocarpus[[7]]
pp_mastocarpus_NW_zi <- pp_mastocarpus[[8]]
pp_mastocarpus_SE_zi <- pp_mastocarpus[[9]]
pp_mastocarpus_SW_zi <- pp_mastocarpus[[10]]

pp_mytilus <- extract(mytilus, model_type = "I_zi")
pp_mytilus_BC <- pp_mytilus[[1]]
pp_mytilus_NE <- pp_mytilus[[2]]
pp_mytilus_NW <- pp_mytilus[[3]]
pp_mytilus_SE <- pp_mytilus[[4]]
pp_mytilus_SW <- pp_mytilus[[5]]
pp_mytilus_BC_zi <- pp_mytilus[[6]]
pp_mytilus_NE_zi <- pp_mytilus[[7]]
pp_mytilus_NW_zi <- pp_mytilus[[8]]
pp_mytilus_SE_zi <- pp_mytilus[[9]]
pp_mytilus_SW_zi <- pp_mytilus[[10]]

pp_nucella <- extract(nucella, model_type = "I_zi")  
pp_nucella_BC <- pp_nucella[[1]]
pp_nucella_NE <- pp_nucella[[2]]
pp_nucella_NW <- pp_nucella[[3]]
pp_nucella_SE <- pp_nucella[[4]]
pp_nucella_SW <- pp_nucella[[5]]
pp_nucella_BC_zi <- pp_nucella[[6]]
pp_nucella_NE_zi <- pp_nucella[[7]]
pp_nucella_NW_zi <- pp_nucella[[8]]
pp_nucella_SE_zi <- pp_nucella[[9]]
pp_nucella_SW_zi <- pp_nucella[[10]]

pp_semibalanus <- extract(semibalanus, model_type = "I_zi")
pp_semibalanus_BC <- pp_semibalanus[[1]]
pp_semibalanus_NE <- pp_semibalanus[[2]]
pp_semibalanus_NW <- pp_semibalanus[[3]]
pp_semibalanus_SE <- pp_semibalanus[[4]]
pp_semibalanus_SW <- pp_semibalanus[[5]]
pp_semibalanus_BC_zi <- pp_semibalanus[[6]]
pp_semibalanus_NE_zi <- pp_semibalanus[[7]]
pp_semibalanus_NW_zi <- pp_semibalanus[[8]]
pp_semibalanus_SE_zi <- pp_semibalanus[[9]]
pp_semibalanus_SW_zi <- pp_semibalanus[[10]]


pp_amphipoda <- extract(amphipoda, model_type = "I_zi")
pp_amphipoda_BC <- pp_amphipoda[[1]]
pp_amphipoda_NE <- pp_amphipoda[[2]]
pp_amphipoda_NW <- pp_amphipoda[[3]]
pp_amphipoda_SE <- pp_amphipoda[[4]]
pp_amphipoda_SW <- pp_amphipoda[[5]]
pp_amphipoda_BC_zi <- pp_amphipoda[[6]]
pp_amphipoda_NE_zi <- pp_amphipoda[[7]]
pp_amphipoda_NW_zi <- pp_amphipoda[[8]]
pp_amphipoda_SE_zi <- pp_amphipoda[[9]]
pp_amphipoda_SW_zi <- pp_amphipoda[[10]]

pp_isopoda <- extract(isopoda, model_type = "I_zi")
pp_isopoda_BC <- pp_isopoda[[1]]
pp_isopoda_NE <- pp_isopoda[[2]]
pp_isopoda_NW <- pp_isopoda[[3]]
pp_isopoda_SE <- pp_isopoda[[4]]
pp_isopoda_SW <- pp_isopoda[[5]]
pp_isopoda_BC_zi <- pp_isopoda[[6]]
pp_isopoda_NE_zi <- pp_isopoda[[7]]
pp_isopoda_NW_zi <- pp_isopoda[[8]]
pp_isopoda_SE_zi <- pp_isopoda[[9]]
pp_isopoda_SW_zi <- pp_isopoda[[10]]

## Function to calculate derivative =============================================
source("Scripts/first_deriv_function.R")

# Calculate second derivative =================================================

deriv1_cover_MDS1_global <- Deriv(global_cover_MDS1)
deriv1_cover_MDS2_global <- Deriv(global_cover_MDS2)
deriv1_cover_MDS3_global <- Deriv(global_cover_MDS3)

deriv1_count_MDS1_global <- Deriv(global_count_MDS1)
deriv1_count_MDS2_global <- Deriv(global_count_MDS2)
deriv1_count_MDS3_global <- Deriv(global_count_MDS3)


deriv1_count_MDS1_NE <- Deriv(pp_site_count_MDS1_NE)
deriv1_count_MDS1_SE <- Deriv(pp_site_count_MDS1_SE)
deriv1_count_MDS1_NW <- Deriv(pp_site_count_MDS1_NW)
deriv1_count_MDS1_SW <- Deriv(pp_site_count_MDS1_SW)
deriv1_count_MDS1_BC <- Deriv(pp_site_count_MDS1_BC)

deriv1_count_MDS2_NE <- Deriv(pp_site_count_MDS2_NE)
deriv1_count_MDS2_SE <- Deriv(pp_site_count_MDS2_SE)
deriv1_count_MDS2_NW <- Deriv(pp_site_count_MDS2_NW)
deriv1_count_MDS2_SW <- Deriv(pp_site_count_MDS2_SW)
deriv1_count_MDS2_BC <- Deriv(pp_site_count_MDS2_BC)

deriv1_count_MDS3_NE <- Deriv(pp_site_count_MDS3_NE)
deriv1_count_MDS3_SE <- Deriv(pp_site_count_MDS3_SE)
deriv1_count_MDS3_NW <- Deriv(pp_site_count_MDS3_NW)
deriv1_count_MDS3_SW <- Deriv(pp_site_count_MDS3_SW)
deriv1_count_MDS3_BC <- Deriv(pp_site_count_MDS3_BC)

deriv1_cover_MDS1_NE <- Deriv(pp_site_cover_MDS1_NE)
deriv1_cover_MDS1_SE <- Deriv(pp_site_cover_MDS1_SE)
deriv1_cover_MDS1_NW <- Deriv(pp_site_cover_MDS1_NW)
deriv1_cover_MDS1_SW <- Deriv(pp_site_cover_MDS1_SW)
deriv1_cover_MDS1_BC <- Deriv(pp_site_cover_MDS1_BC)

deriv1_cover_MDS2_NE <- Deriv(pp_site_cover_MDS2_NE)
deriv1_cover_MDS2_SE <- Deriv(pp_site_cover_MDS2_SE)
deriv1_cover_MDS2_NW <- Deriv(pp_site_cover_MDS2_NW)
deriv1_cover_MDS2_SW <- Deriv(pp_site_cover_MDS2_SW)
deriv1_cover_MDS2_BC <- Deriv(pp_site_cover_MDS2_BC)

deriv1_cover_MDS3_NE <- Deriv(pp_site_cover_MDS3_NE)
deriv1_cover_MDS3_SE <- Deriv(pp_site_cover_MDS3_SE)
deriv1_cover_MDS3_NW <- Deriv(pp_site_cover_MDS3_NW)
deriv1_cover_MDS3_SW <- Deriv(pp_site_cover_MDS3_SW)
deriv1_cover_MDS3_BC <- Deriv(pp_site_cover_MDS3_BC)


deriv1_ascophylum_BC <- Deriv(pp_ascophylum_BC)
deriv1_ascophylum_NE <- Deriv(pp_ascophylum_NE)
deriv1_ascophylum_NW <- Deriv(pp_ascophylum_NW)
deriv1_ascophylum_SE <- Deriv(pp_ascophylum_SE)
deriv1_ascophylum_SW <- Deriv(pp_ascophylum_SW)
deriv1_ascophylum_BC_zi <- Deriv(pp_ascophylum_BC_zi)
deriv1_ascophylum_NE_zi <- Deriv(pp_ascophylum_NE_zi)
deriv1_ascophylum_NW_zi <- Deriv(pp_ascophylum_NW_zi)
deriv1_ascophylum_SE_zi <- Deriv(pp_ascophylum_SE_zi)
deriv1_ascophylum_SW_zi <- Deriv(pp_ascophylum_SW_zi)


deriv1_chondrus_BC <- Deriv(pp_chondrus_BC)
deriv1_chondrus_NE <- Deriv(pp_chondrus_NE)
deriv1_chondrus_NW <- Deriv(pp_chondrus_NW)
deriv1_chondrus_SE <- Deriv(pp_chondrus_SE)
deriv1_chondrus_SW <- Deriv(pp_chondrus_SW)
deriv1_chondrus_BC_zi <- Deriv(pp_chondrus_BC_zi)
deriv1_chondrus_NE_zi <- Deriv(pp_chondrus_NE_zi)
deriv1_chondrus_NW_zi <- Deriv(pp_chondrus_NW_zi)
deriv1_chondrus_SE_zi <- Deriv(pp_chondrus_SE_zi)
deriv1_chondrus_SW_zi <- Deriv(pp_chondrus_SW_zi)

deriv1_cmaenas_BC <- Deriv(pp_cmaenas_BC) 
deriv1_cmaenas_NE <- Deriv(pp_cmaenas_NE)
deriv1_cmaenas_NW <- Deriv(pp_cmaenas_NW)
deriv1_cmaenas_SE <- Deriv(pp_cmaenas_SE)
deriv1_cmaenas_SW <- Deriv(pp_cmaenas_SW)
deriv1_cmaenas_BC_zi <- Deriv(pp_cmaenas_BC_zi) 
deriv1_cmaenas_NE_zi <- Deriv(pp_cmaenas_NE_zi)
deriv1_cmaenas_NW_zi <- Deriv(pp_cmaenas_NW_zi)
deriv1_cmaenas_SE_zi <- Deriv(pp_cmaenas_SE_zi)
deriv1_cmaenas_SW_zi <- Deriv(pp_cmaenas_SW_zi)

deriv1_fucus_BC <- Deriv(pp_fucus_BC)
deriv1_fucus_NE <- Deriv(pp_fucus_NE)
deriv1_fucus_NW <- Deriv(pp_fucus_NW)
deriv1_fucus_SE <- Deriv(pp_fucus_SE)
deriv1_fucus_SW <- Deriv(pp_fucus_SW)
deriv1_fucus_BC_zi <- Deriv(pp_fucus_BC_zi)
deriv1_fucus_NE_zi <- Deriv(pp_fucus_NE_zi)
deriv1_fucus_NW_zi <- Deriv(pp_fucus_NW_zi)
deriv1_fucus_SE_zi <- Deriv(pp_fucus_SE_zi)
deriv1_fucus_SW_zi <- Deriv(pp_fucus_SW_zi)

deriv1_littorina_BC <- Deriv(pp_littorina_BC)
deriv1_littorina_NE <- Deriv(pp_littorina_NE)
deriv1_littorina_NW <- Deriv(pp_littorina_NW)
deriv1_littorina_SE <- Deriv(pp_littorina_SE)
deriv1_littorina_SW <- Deriv(pp_littorina_SW)
deriv1_littorina_BC_zi <- Deriv(pp_littorina_BC_zi)
deriv1_littorina_NE_zi <- Deriv(pp_littorina_NE_zi)
deriv1_littorina_NW_zi <- Deriv(pp_littorina_NW_zi)
deriv1_littorina_SE_zi <- Deriv(pp_littorina_SE_zi)
deriv1_littorina_SW_zi <- Deriv(pp_littorina_SW_zi)

deriv1_mastocarpus_BC <- Deriv(pp_mastocarpus_BC)
deriv1_mastocarpus_NE <- Deriv(pp_mastocarpus_NE)
deriv1_mastocarpus_NW <- Deriv(pp_mastocarpus_NW)
deriv1_mastocarpus_SE <- Deriv(pp_mastocarpus_SE)
deriv1_mastocarpus_SW <- Deriv(pp_mastocarpus_SW)
deriv1_mastocarpus_BC_zi <- Deriv(pp_mastocarpus_BC_zi)
deriv1_mastocarpus_NE_zi <- Deriv(pp_mastocarpus_NE_zi)
deriv1_mastocarpus_NW_zi <- Deriv(pp_mastocarpus_NW_zi)
deriv1_mastocarpus_SE_zi <- Deriv(pp_mastocarpus_SE_zi)
deriv1_mastocarpus_SW_zi <- Deriv(pp_mastocarpus_SW_zi)

deriv1_mytilus_BC <- Deriv(pp_mytilus_BC)
deriv1_mytilus_NE <- Deriv(pp_mytilus_NE)
deriv1_mytilus_NW <- Deriv(pp_mytilus_NW)
deriv1_mytilus_SE <- Deriv(pp_mytilus_SE)
deriv1_mytilus_SW <- Deriv(pp_mytilus_SW)
deriv1_mytilus_BC_zi <- Deriv(pp_mytilus_BC_zi)
deriv1_mytilus_NE_zi <- Deriv(pp_mytilus_NE_zi)
deriv1_mytilus_NW_zi <- Deriv(pp_mytilus_NW_zi)
deriv1_mytilus_SE_zi <- Deriv(pp_mytilus_SE_zi)
deriv1_mytilus_SW_zi <- Deriv(pp_mytilus_SW_zi)

deriv1_nucella_BC <- Deriv(pp_nucella_BC)
deriv1_nucella_NE <- Deriv(pp_nucella_NE)
deriv1_nucella_NW <- Deriv(pp_nucella_NW)
deriv1_nucella_SE <- Deriv(pp_nucella_SE)
deriv1_nucella_SW <- Deriv(pp_nucella_SW)
deriv1_nucella_BC_zi <- Deriv(pp_nucella_BC_zi)
deriv1_nucella_NE_zi <- Deriv(pp_nucella_NE_zi)
deriv1_nucella_NW_zi <- Deriv(pp_nucella_NW_zi)
deriv1_nucella_SE_zi <- Deriv(pp_nucella_SE_zi)
deriv1_nucella_SW_zi <- Deriv(pp_nucella_SW_zi)

deriv1_semibalanus_BC <- Deriv(pp_semibalanus_BC)
deriv1_semibalanus_NE <- Deriv(pp_semibalanus_NE)
deriv1_semibalanus_NW <- Deriv(pp_semibalanus_NW)
deriv1_semibalanus_SE <- Deriv(pp_semibalanus_SE)
deriv1_semibalanus_SW <- Deriv(pp_semibalanus_SW)
deriv1_semibalanus_BC_zi <- Deriv(pp_semibalanus_BC_zi)
deriv1_semibalanus_NE_zi <- Deriv(pp_semibalanus_NE_zi)
deriv1_semibalanus_NW_zi <- Deriv(pp_semibalanus_NW_zi)
deriv1_semibalanus_SE_zi <- Deriv(pp_semibalanus_SE_zi)
deriv1_semibalanus_SW_zi <- Deriv(pp_semibalanus_SW_zi)

deriv1_amphipoda_BC <- Deriv(pp_amphipoda_BC)
deriv1_amphipoda_NE <- Deriv(pp_amphipoda_NE)
deriv1_amphipoda_NW <- Deriv(pp_amphipoda_NW)
deriv1_amphipoda_SE <- Deriv(pp_amphipoda_SE)
deriv1_amphipoda_SW <- Deriv(pp_amphipoda_SW)
deriv1_amphipoda_BC_zi <- Deriv(pp_amphipoda_BC_zi)
deriv1_amphipoda_NE_zi <- Deriv(pp_amphipoda_NE_zi)
deriv1_amphipoda_NW_zi <- Deriv(pp_amphipoda_NW_zi)
deriv1_amphipoda_SE_zi <- Deriv(pp_amphipoda_SE_zi)
deriv1_amphipoda_SW_zi <- Deriv(pp_amphipoda_SW_zi)


deriv1_isopoda_BC <- Deriv(pp_isopoda_BC)
deriv1_isopoda_NE <- Deriv(pp_isopoda_NE)
deriv1_isopoda_NW <- Deriv(pp_isopoda_NW)
deriv1_isopoda_SE <- Deriv(pp_isopoda_SE)
deriv1_isopoda_SW <- Deriv(pp_isopoda_SW)
deriv1_isopoda_BC_zi <- Deriv(pp_isopoda_BC_zi)
deriv1_isopoda_NE_zi <- Deriv(pp_isopoda_NE_zi)
deriv1_isopoda_NW_zi <- Deriv(pp_isopoda_NW_zi)
deriv1_isopoda_SE_zi <- Deriv(pp_isopoda_SE_zi)
deriv1_isopoda_SW_zi <- Deriv(pp_isopoda_SW_zi)


##### Make table with first deriv values for each specie/axis and regime shift
#============================================================================

values_deriv1 <- function(data = deriv1_cover_MDS1, rs1 = c(1987,1992.7),
                          rs2 = c(1999.4,2002), rs3 = c(2009.5,2011.7),
                          name = "Sessile MDS1"){
  # create a df for future output
  
  df <- data.frame(var = c(name,name,name),
                   regime_shift = c(1,2,3),
                   upper = c(NA,NA,NA),
                   lower = c(NA,NA,NA))
  
  # make 1 df for data of each regime shift
  data1 <- filter(data, rs1[1] <= data$year)
  data1 <- filter(data1, data1$year <= rs1[2])
  
  data2 <- filter(data, rs2[1] <= data$year)
  data2 <- filter(data2, data2$year <= rs2[2])
  
  data3 <- filter(data, rs3[1] <= data$year)
  data3 <- filter(data3, data3$year <= rs3[2])
  
  # make a column with upper*lower so the max is where both are + or -
  
  data1$vector = data1$upper_CI * data1$lower_CI
  data2$vector = data2$upper_CI * data2$lower_CI
  data3$vector = data3$upper_CI * data3$lower_CI
  
  # finding where vector is max
  
  ind1 <- which.max(data1$vector)
  ind2 <- which.max(data2$vector)
  ind3 <- which.max(data3$vector)
  
  if (data1$vector[ind1] > 0){
    df$upper[1] <- data1$upper_CI[ind1]
    df$lower[1] <- data1$lower_CI[ind1]
  } else{
    ind <- which.max(abs(data1$y_deriv_estim))
    df$upper[1] <- data1$upper_CI[ind]
    df$lower[1] <- data1$lower_CI[ind]
  }
  
  if (data2$vector[ind2] > 0){
    df$upper[2] <- data2$upper_CI[ind2]
    df$lower[2] <- data2$lower_CI[ind2]
  } else{
    ind <- which.max(abs(data2$y_deriv_estim))
    df$upper[2] <- data2$upper_CI[ind]
    df$lower[2] <- data2$lower_CI[ind]
  }
  
  if (data3$vector[ind3] > 0){
    df$upper[3] <- data3$upper_CI[ind3]
    df$lower[3] <- data3$lower_CI[ind3]
  } else{
    ind <- which.max(abs(data3$y_deriv_estim))
    df$upper[3] <- data3$upper_CI[ind]
    df$lower[3] <- data3$lower_CI[ind]
  }
  return(df)  
  
}

#apply the function for all second deriv and merge in 1 DF 

vm1 <- values_deriv1(data=deriv1_cover_MDS1_global, name= "Sessile MDS1 global")
vm2 <- values_deriv1(data=deriv1_cover_MDS2_global, name= "Sessile MDS2 global")
vm3 <- values_deriv1(data=deriv1_cover_MDS3_global, name= "Sessile MDS2 global")

vm4 <- values_deriv1(data=deriv1_count_MDS1_global, name= "Mobile MDS1 global")
vm5 <- values_deriv1(data=deriv1_count_MDS2_global, name= "Mobile MDS2 global")
vm6 <- values_deriv1(data=deriv1_count_MDS3_global, name= "Mobile MDS2 global")

vm7 <- values_deriv1(data=deriv1_count_MDS1_BC, name = "Mobile MDS1 BC")
vm8 <- values_deriv1(data=deriv1_count_MDS1_NE, name = "Mobile MDS1 NE")
vm9 <- values_deriv1(data=deriv1_count_MDS1_NW, name = "Mobile MDS1 NW")
vm10 <- values_deriv1(data=deriv1_count_MDS1_SE, name = "Mobile MDS1 SE")
vm11 <- values_deriv1(data=deriv1_count_MDS1_SW, name = "Mobile MDS1 SW")

vm7 <- values_deriv1(data=deriv1_count_MDS2_BC, name = "Mobile MDS2 BC")
vm8 <- values_deriv1(data=deriv1_count_MDS2_NE, name = "Mobile MDS2 NE")
vm9 <- values_deriv1(data=deriv1_count_MDS2_NW, name = "Mobile MDS2 NW")
vm10 <- values_deriv1(data=deriv1_count_MDS2_SE, name = "Mobile MDS2 SE")
vm11 <- values_deriv1(data=deriv1_count_MDS2_SW, name = "Mobile MDS2 SW")

vm12 <- values_deriv1(data=deriv1_count_MDS3_BC, name = "Mobile MDS3 BC")
vm13 <- values_deriv1(data=deriv1_count_MDS3_NE, name = "Mobile MDS3 NE")
vm14 <- values_deriv1(data=deriv1_count_MDS3_NW, name = "Mobile MDS3 NW")
vm15 <- values_deriv1(data=deriv1_count_MDS3_SE, name = "Mobile MDS3 SE")
vm16 <- values_deriv1(data=deriv1_count_MDS3_SW, name = "Mobile MDS3 SW")

vm17 <- values_deriv1(data=deriv1_cover_MDS1_BC, name = "Sessile MDS1 BC")
vm18 <- values_deriv1(data=deriv1_cover_MDS1_NE, name = "Sessile MDS1 NE")
vm19 <- values_deriv1(data=deriv1_cover_MDS1_NW, name = "Sessile MDS1 NW")
vm20 <- values_deriv1(data=deriv1_cover_MDS1_SE, name = "Sessile MDS1 SE")
vm21 <- values_deriv1(data=deriv1_cover_MDS1_SW, name = "Sessile MDS1 SW")

vm22 <- values_deriv1(data=deriv1_cover_MDS2_BC, name = "Sessile MDS2 BC")
vm23 <- values_deriv1(data=deriv1_cover_MDS2_NE, name = "Sessile MDS2 NE")
vm24 <- values_deriv1(data=deriv1_cover_MDS2_NW, name = "Sessile MDS2 NW")
vm25 <- values_deriv1(data=deriv1_cover_MDS2_SE, name = "Sessile MDS2 SE")
vm26 <- values_deriv1(data=deriv1_cover_MDS2_SW, name = "Sessile MDS2 SW")

vm27 <- values_deriv1(data=deriv1_cover_MDS3_BC, name = "Sessile MDS3 BC")
vm28 <- values_deriv1(data=deriv1_cover_MDS3_NE, name = "Sessile MDS3 NE")
vm29 <- values_deriv1(data=deriv1_cover_MDS3_NW, name = "Sessile MDS3 NW")
vm30 <- values_deriv1(data=deriv1_cover_MDS3_SE, name = "Sessile MDS3 SE")
vm31 <- values_deriv1(data=deriv1_cover_MDS3_SW, name = "Sessile MDS3 SW")


v11 <- values_deriv1(data = deriv1_ascophylum_BC, name = "Ascophyllum BC")
v12 <- values_deriv1(data = deriv1_ascophylum_NE, name = "Ascophyllum NE")
v13 <- values_deriv1(data = deriv1_ascophylum_NW, name = "Ascophyllum NW")
v14 <- values_deriv1(data = deriv1_ascophylum_SE, name = "Ascophyllum SE")
v15 <- values_deriv1(data = deriv1_ascophylum_SW, name = "Ascophyllum SW")

v16 <- values_deriv1(data = deriv1_chondrus_BC, name = "Chondrus BC")
v17 <- values_deriv1(data = deriv1_chondrus_NE, name = "Chondrus NE")
v18 <- values_deriv1(data = deriv1_chondrus_NW, name = "Chondrus NW")
v19 <- values_deriv1(data = deriv1_chondrus_SE, name = "Chondrus SE")
v20 <- values_deriv1(data = deriv1_chondrus_SW, name = "Chondrus SW")

v21 <- values_deriv1(data = deriv1_cmaenas_BC, name = "C. maenas BC")
v22 <- values_deriv1(data = deriv1_cmaenas_NE, name = "C. maenas NE")
v23 <- values_deriv1(data = deriv1_cmaenas_NW, name = "C. maenas NW")
v24 <- values_deriv1(data = deriv1_cmaenas_SE, name = "C. maenas SE")
v25 <- values_deriv1(data = deriv1_cmaenas_SW, name = "C. maenas SW")

v26 <- values_deriv1(data = deriv1_fucus_BC, name = "Fucus BC")
v27 <- values_deriv1(data = deriv1_fucus_NE, name = "Fucus NE")
v28 <- values_deriv1(data = deriv1_fucus_NW, name = "Fucus NW")
v29 <- values_deriv1(data = deriv1_fucus_SE, name = "Fucus SE")
v30 <- values_deriv1(data = deriv1_fucus_SW, name = "Fucus SW")

v31 <- values_deriv1(data = deriv1_littorina_BC, name = "Littorina BC")
v32 <- values_deriv1(data = deriv1_littorina_NE, name = "Littorina NE")
v33 <- values_deriv1(data = deriv1_littorina_NW, name = "Littorina NW")
v34 <- values_deriv1(data = deriv1_littorina_SE, name = "Littorina SE")
v35 <- values_deriv1(data = deriv1_littorina_SW, name = "Littorina SW")

v36 <- values_deriv1(data = deriv1_mastocarpus_BC, name = "Mastocarpus BC")
v37 <- values_deriv1(data = deriv1_mastocarpus_NE, name = "Mastocarpus NE")
v38 <- values_deriv1(data = deriv1_mastocarpus_NW, name = "Mastocarpus NW")
v39 <- values_deriv1(data = deriv1_mastocarpus_SE, name = "Mastocarpus SE")
v40 <- values_deriv1(data = deriv1_mastocarpus_SW, name = "Mastocarpus SW")

v41 <- values_deriv1(data = deriv1_mytilus_BC, name = "mytilus BC")
v42 <- values_deriv1(data = deriv1_mytilus_NE, name = "Mytilus NE")
v43 <- values_deriv1(data = deriv1_mytilus_NW, name = "Mytilus NW")
v44 <- values_deriv1(data = deriv1_mytilus_SE, name = "Mytilus SE")
v45 <- values_deriv1(data = deriv1_mytilus_SW, name = "Mytilus SW")

v46 <- values_deriv1(data = deriv1_nucella_BC, name = "Nucella BC")
v47 <- values_deriv1(data = deriv1_nucella_NE, name = "Nucella NE")
v48 <- values_deriv1(data = deriv1_nucella_NW, name = "Nucella NW")
v49 <- values_deriv1(data = deriv1_nucella_SE, name = "Nucella SE")
v50 <- values_deriv1(data = deriv1_nucella_SW, name = "Nucella SW")

v51 <- values_deriv1(data = deriv1_semibalanus_BC, name = "Semibalanus BC")
v52 <- values_deriv1(data = deriv1_semibalanus_NE, name = "Semibalanus NE")
v53 <- values_deriv1(data = deriv1_semibalanus_NW, name = "Semibalanus NW")
v54 <- values_deriv1(data = deriv1_semibalanus_SE, name = "Semibalanus SE")
v55 <- values_deriv1(data = deriv1_semibalanus_SW, name = "Semibalanus SW")

v56 <- values_deriv1(data = deriv1_amphipoda_BC, name = "Amphipoda BC")
v57 <- values_deriv1(data = deriv1_amphipoda_NE, name = "Amphipoda NE")
v58 <- values_deriv1(data = deriv1_amphipoda_NW, name = "Amphipoda NW")
v59 <- values_deriv1(data = deriv1_amphipoda_SE, name = "Amphipoda SE")
v60 <- values_deriv1(data = deriv1_amphipoda_SW, name = "Amphipoda SW")

v61 <- values_deriv1(data = deriv1_isopoda_BC, name = "Isopoda BC")
v62 <- values_deriv1(data = deriv1_isopoda_NE, name = "Isopoda NE")
v63 <- values_deriv1(data = deriv1_isopoda_NW, name = "Isopoda NW")
v64 <- values_deriv1(data = deriv1_isopoda_SE, name = "Isopoda SE")
v65 <- values_deriv1(data = deriv1_isopoda_SW, name = "Isopoda SW")


table1 <- rbind(vm1,vm2,vm3,vm4,vm6,vm5,vm7,vm8,vm9,vm10,
                vm11,vm12,vm13,vm14,vm15,vm16,vm17,vm18,vm19,
                vm20,vm21,vm22,vm23,vm24,vm25,vm26,vm27,vm28,
                vm29,vm30,vm31,
                v11,v12,v13,v14,v15,v16,v17,v18,v19,
                v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31,v32,v33,v34,v35,v36,
                v37,v38,v39,v40,v41,v42,v43,v44,v45,v46,v47,v48,v49,v50,v51,v52,v53
                ,v54,v55,v56,v57,v58,v59,v60,v61,v62,v63,v64,v65)


### Same for the zi and merge the two

vm1$upper <- NA # replace per NA for models where there is no zi
vm1$lower <- NA
vm2$upper <- NA
vm2$lower <- NA
vm3$upper <- NA
vm3$lower <- NA
vm4$upper <- NA
vm4$lower <- NA
vm5$upper <- NA
vm5$lower <- NA
vm7$upper <- NA
vm7$lower <- NA
vm8$upper <- NA
vm8$lower <- NA
vm9$upper <- NA
vm9$lower <- NA
vm10$upper <- NA
vm10$lower <- NA
vm11$upper <- NA # replace per NA for models where there is no zi
vm11$lower <- NA
vm12$upper <- NA
vm12$lower <- NA
vm13$upper <- NA
vm13$lower <- NA
vm14$upper <- NA
vm14$lower <- NA
vm15$upper <- NA
vm15$lower <- NA
vm17$upper <- NA
vm17$lower <- NA
vm16$lower <-NA
vm16$upper <-NA
vm18$upper <- NA
vm18$lower <- NA
vm19$upper <- NA
vm19$lower <- NA
vm20$upper <- NA
vm20$lower <- NA
vm21$upper <- NA # replace per NA for models where there is no zi
vm21$lower <- NA
vm22$upper <- NA
vm22$lower <- NA
vm23$upper <- NA
vm23$lower <- NA
vm24$upper <- NA
vm24$lower <- NA
vm25$upper <- NA
vm25$lower <- NA
vm27$upper <- NA
vm27$lower <- NA
vm26$upper <- NA
vm26$lower <- NA
vm28$upper <- NA
vm28$lower <- NA
vm29$upper <- NA
vm29$lower <- NA
vm30$upper <- NA
vm30$lower <- NA
vm31$upper <- NA
vm31$lower <- NA

v11 <- values_deriv1(data = deriv1_ascophylum_BC_zi, name = "Ascophyllum BC")
v12 <- values_deriv1(data = deriv1_ascophylum_NE_zi, name = "Ascophyllum NE")
v13 <- values_deriv1(data = deriv1_ascophylum_NW_zi, name = "Ascophyllum NW")
v14 <- values_deriv1(data = deriv1_ascophylum_SE_zi, name = "Ascophyllum SE")
v15 <- values_deriv1(data = deriv1_ascophylum_SW, name = "Ascophyllum SW")

v16 <- values_deriv1(data = deriv1_chondrus_BC_zi, name = "Chondrus BC")
v17 <- values_deriv1(data = deriv1_chondrus_NE_zi, name = "Chondrus NE")
v18 <- values_deriv1(data = deriv1_chondrus_NW_zi, name = "Chondrus NW")
v19 <- values_deriv1(data = deriv1_chondrus_SE_zi, name = "Chondrus SE")
v20 <- values_deriv1(data = deriv1_chondrus_SW_zi, name = "Chondrus SW")

v21 <- values_deriv1(data = deriv1_cmaenas_BC_zi, name = "C. maenas BC")
v22 <- values_deriv1(data = deriv1_cmaenas_NE_zi, name = "C. maenas NE")
v23 <- values_deriv1(data = deriv1_cmaenas_NW_zi, name = "C. maenas NW")
v24 <- values_deriv1(data = deriv1_cmaenas_SE_zi, name = "C. maenas SE")
v25 <- values_deriv1(data = deriv1_cmaenas_SW_zi, name = "C. maenas SW")

v26 <- values_deriv1(data = deriv1_fucus_BC_zi, name = "Fucus BC")
v27 <- values_deriv1(data = deriv1_fucus_NE_zi, name = "Fucus NE")
v28 <- values_deriv1(data = deriv1_fucus_NW_zi, name = "Fucus NW")
v29 <- values_deriv1(data = deriv1_fucus_SE_zi, name = "Fucus SE")
v30 <- values_deriv1(data = deriv1_fucus_SW_zi, name = "Fucus SW")

v31 <- values_deriv1(data = deriv1_littorina_BC_zi, name = "Littorina BC")
v32 <- values_deriv1(data = deriv1_littorina_NE_zi, name = "Littorina NE")
v33 <- values_deriv1(data = deriv1_littorina_NW_zi, name = "Littorina NW")
v34 <- values_deriv1(data = deriv1_littorina_SE_zi, name = "Littorina SE")
v35 <- values_deriv1(data = deriv1_littorina_SW_zi, name = "Littorina SW")

v36 <- values_deriv1(data = deriv1_mastocarpus_BC_zi, name = "Mastocarpus BC")
v37 <- values_deriv1(data = deriv1_mastocarpus_NE_zi, name = "Mastocarpus NE")
v38 <- values_deriv1(data = deriv1_mastocarpus_NW_zi, name = "Mastocarpus NW")
v39 <- values_deriv1(data = deriv1_mastocarpus_SE_zi, name = "Mastocarpus SE")
v40 <- values_deriv1(data = deriv1_mastocarpus_SW_zi, name = "Mastocarpus SW")

v41 <- values_deriv1(data = deriv1_mytilus_BC_zi, name = "mytilus BC")
v42 <- values_deriv1(data = deriv1_mytilus_NE_zi, name = "Mytilus NE")
v43 <- values_deriv1(data = deriv1_mytilus_NW_zi, name = "Mytilus NW")
v44 <- values_deriv1(data = deriv1_mytilus_SE_zi, name = "Mytilus SE")
v45 <- values_deriv1(data = deriv1_mytilus_SW_zi, name = "Mytilus SW")

v46 <- values_deriv1(data = deriv1_nucella_BC_zi, name = "Nucella BC")
v47 <- values_deriv1(data = deriv1_nucella_NE_zi, name = "Nucella NE")
v48 <- values_deriv1(data = deriv1_nucella_NW_zi, name = "Nucella NW")
v49 <- values_deriv1(data = deriv1_nucella_SE_zi, name = "Nucella SE")
v50 <- values_deriv1(data = deriv1_nucella_SW_zi, name = "Nucella SW")

v51 <- values_deriv1(data = deriv1_semibalanus_BC_zi, name = "Semibalanus BC")
v52 <- values_deriv1(data = deriv1_semibalanus_NE_zi, name = "Semibalanus NE")
v53 <- values_deriv1(data = deriv1_semibalanus_NW_zi, name = "Semibalanus NW")
v54 <- values_deriv1(data = deriv1_semibalanus_SE_zi, name = "Semibalanus SE")
v55 <- values_deriv1(data = deriv1_semibalanus_SW_zi, name = "Semibalanus SW")


v56 <- values_deriv1(data = deriv1_amphipoda_BC_zi, name = "Amphipoda BC")
v57 <- values_deriv1(data = deriv1_amphipoda_NE_zi, name = "Amphipoda NE")
v58 <- values_deriv1(data = deriv1_amphipoda_NW_zi, name = "Amphipoda NW")
v59 <- values_deriv1(data = deriv1_amphipoda_SE_zi, name = "Amphipoda SE")
v60 <- values_deriv1(data = deriv1_amphipoda_SW_zi, name = "Amphipoda SW")

v61 <- values_deriv1(data = deriv1_isopoda_BC_zi, name = "Isopoda BC")
v62 <- values_deriv1(data = deriv1_isopoda_NE_zi, name = "Isopoda NE")
v63 <- values_deriv1(data = deriv1_isopoda_NW_zi, name = "Isopoda NW")
v64 <- values_deriv1(data = deriv1_isopoda_SE_zi, name = "Isopoda SE")
v65 <- values_deriv1(data = deriv1_isopoda_SW_zi, name = "Isopoda SW")



table2 <- rbind(vm1,vm2,vm3,vm4,vm6,vm5,vm7,vm8,vm9,vm10,
                vm11,vm12,vm13,vm14,vm15,vm16,vm17,vm18,vm19,
                vm20,vm21,vm22,vm23,vm24,vm25,vm26,vm27,vm28,
                vm29,vm30,vm31,
                v11,v12,v13,v14,v15,v16,v17,v18,v19,
                v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31,v32,v33,v34,v35,v36,
                v37,v38,v39,v40,v41,v42,v43,v44,v45,v46,v47,v48,v49,v50,v51,v52,v53
                ,v54,v55,v56,v57,v58,v59,v60,v61,v62,v63,v64,v65)


colnames(table2) <- c("var","regime_shift","upper_zi","lower_zi")

table1 <- unique(merge(table1,table2, by=c("var","regime_shift")))


# # add R2 values 
# 
# table1$R2 <- NA
# 
# table1$R2[1*3-2] <- mean(cover$criteria$bayes_R2[,1])
# table1$R2[2*3-2] <- mean(cover$criteria$bayes_R2[,2])
# 
# table1$R2[3*3-2] <- mean(count_4$criteria$bayes_R2[,1])
# table1$R2[4*3-2] <- mean(count_4$criteria$bayes_R2[,2])
# 
# table1$R2[11*3-2] <- mean(ascophyllum$criteria$bayes_R2)
# 
# table1$R2[12*3-2] <- mean(chondrus$criteria$bayes_R2)
# 
# table1$R2[18*3-2] <- mean(cmenas$criteria$bayes_R2)
# 
# table1$R2[24*3-2] <- mean(fucus$criteria$bayes_R2)
# 
# table1$R2[25*3-2] <- mean(littorina$criteria$bayes_R2)
# 
# table1$R2[31*3-2] <- mean(mastocarpus$criteria$bayes_R2)
# 
# table1$R2[37*3-2] <- mean(mytilus$criteria$bayes_R2)
# 
# table1$R2[43*3-2] <- mean(nucella$criteria$bayes_R2)
# 
# table1$R2[49*3-2] <- mean(semibalanus$criteria$bayes_R2)

write.csv(table1,file = "time_series_outputs/values_deriv1.csv")



### Function to plot derivative ######################################
source("Scripts/plot_deriv.R")

cover_MDS1_global <- plot_deriv(deriv1_cover_MDS1_global,"deriv1(cover MDS1 global)")
cover_MDS2_global <- plot_deriv(deriv1_cover_MDS2_global,"deriv1(cover MDS2 global)")
cover_MDS3_global <- plot_deriv(deriv1_cover_MDS3_global,"deriv1(cover MDS3 global)")

count_MDS1_global <- plot_deriv(deriv1_count_MDS1_global,"deriv1(count MDS1 global)")
count_MDS2_global <- plot_deriv(deriv1_count_MDS2_global,"deriv1(count MDS2 global)")
count_MDS3_global <- plot_deriv(deriv1_count_MDS3_global,"deriv1(count MDS3 global)")


count_MDS1_BC <- plot_deriv(deriv1_count_MDS1_BC, "deriv1(count MDS1 BC)") # to fix
count_MDS1_SE <- plot_deriv(deriv1_count_MDS1_SE, "deriv1(count MDS1 SE)") 
count_MDS1_NW <- plot_deriv(deriv1_count_MDS1_NW, "deriv1(count MDS1 NW)") 
count_MDS1_SW <- plot_deriv(deriv1_count_MDS1_SW, "deriv1(count MDS1 SW")
count_MDS1_NE <- plot_deriv(deriv1_count_MDS1_NE, "deriv1(count MDS1 NE")

count_MDS2_BC <- plot_deriv(deriv1_count_MDS2_BC, "deriv1(count MDS2 BC)") # to fix
count_MDS2_SE <- plot_deriv(deriv1_count_MDS2_SE, "deriv1(count MDS2 SE)") 
count_MDS2_NW <- plot_deriv(deriv1_count_MDS2_NW, "deriv1(count MDS2 NW)") 
count_MDS2_SW <- plot_deriv(deriv1_count_MDS2_SW, "deriv1(count MDS2 SW)")
count_MDS2_NE <- plot_deriv(deriv1_count_MDS2_NE, "deriv1(count MDS2 NE)")

count_MDS3_BC <- plot_deriv(deriv1_count_MDS3_BC, "deriv1(count MDS3 BC)") # to fix
count_MDS3_SE <- plot_deriv(deriv1_count_MDS3_SE, "deriv1(count MDS3 SE)") 
count_MDS3_NW <- plot_deriv(deriv1_count_MDS3_NW, "deriv1(count MDS3 NW)") 
count_MDS3_SW <- plot_deriv(deriv1_count_MDS3_SW, "deriv1(count MDS3 SW")
count_MDS3_NE <- plot_deriv(deriv1_count_MDS3_NE, "deriv1(count MDS3 NE")

cover_MDS1_BC <- plot_deriv(deriv1_cover_MDS1_BC, "deriv1(cover MDS1 BC)") # to fix
cover_MDS1_SE <- plot_deriv(deriv1_cover_MDS1_SE, "deriv1(cover MDS1 SE)") 
cover_MDS1_NW <- plot_deriv(deriv1_cover_MDS1_NW, "deriv1(cover MDS1 NW)") 
cover_MDS1_SW <- plot_deriv(deriv1_cover_MDS1_SW, "deriv1(cover MDS1 SW")
cover_MDS1_NE <- plot_deriv(deriv1_cover_MDS1_NE, "deriv1(cover MDS1 NE")

cover_MDS2_BC <- plot_deriv(deriv1_cover_MDS2_BC, "deriv1(cover MDS2 BC)") # to fix
cover_MDS2_SE <- plot_deriv(deriv1_cover_MDS2_SE, "deriv1(cover MDS2 SE)") 
cover_MDS2_NW <- plot_deriv(deriv1_cover_MDS2_NW, "deriv1(cover MDS2 NW)") 
cover_MDS2_SW <- plot_deriv(deriv1_cover_MDS2_SW, "deriv1(cover MDS2 SW)")
cover_MDS2_NE <- plot_deriv(deriv1_cover_MDS2_NE, "deriv1(cover MDS2 NE)")

cover_MDS3_BC <- plot_deriv(deriv1_cover_MDS3_BC, "deriv1(cover MDS3 BC)") # to fix
cover_MDS3_SE <- plot_deriv(deriv1_cover_MDS3_SE, "deriv1(cover MDS3 SE)") 
cover_MDS3_NW <- plot_deriv(deriv1_cover_MDS3_NW, "deriv1(cover MDS3 NW)") 
cover_MDS3_SW <- plot_deriv(deriv1_cover_MDS3_SW, "deriv1(cover MDS3 SW)")
cover_MDS3_NE <- plot_deriv(deriv1_cover_MDS3_NE, "deriv1(cover MDS3 NE)")

ascophyllum_BC <- plot_deriv(deriv1_ascophylum_BC, "deriv1(ascophyllum BC)")
ascophyllum_NE <- plot_deriv(deriv1_ascophylum_NE, "deriv1(ascophyllum NE)")
ascophyllum_NW <- plot_deriv(deriv1_ascophylum_BC, "deriv1(ascophyllum NW)")
ascophyllum_SE <- plot_deriv(deriv1_ascophylum_NE, "deriv1(ascophyllum SE)")
ascophyllum_SW <- plot_deriv(deriv1_ascophylum_BC, "deriv1(ascophyllum SW)")

chondrus_BC <- plot_deriv(deriv1_chondrus_BC, "deriv1(chondrus BC)")
#chondrus_MC <- plot_deriv(deriv1_chondrus_MC, "deriv1(chondrus MC)")
chondrus_NE <- plot_deriv(deriv1_chondrus_NE, "deriv1(chondrus NE)")
chondrus_NW <- plot_deriv(deriv1_chondrus_NW, "deriv1(chondrus NW)")
chondrus_SE <- plot_deriv(deriv1_chondrus_SE, "deriv1(chondrus SE)")
chondrus_SW <- plot_deriv(deriv1_chondrus_SW, "deriv1(chondrus SW)")

cmaenas_BC <- plot_deriv(deriv1_cmaenas_BC, "deriv1(C maenas BC)")
#cmaenas_MC <- plot_deriv(deriv1_cmaenas_MC, "deriv1(C maenas MC)")
cmaenas_NE <- plot_deriv(deriv1_cmaenas_NE, "deriv1(C maenas NE)")
cmaenas_NW <- plot_deriv(deriv1_cmaenas_NW, "deriv1(C maenas NW)")
cmaenas_SE <- plot_deriv(deriv1_cmaenas_SE, "deriv1(C maenas SE)")
cmaenas_SW <- plot_deriv(deriv1_cmaenas_SW, "deriv1(C maenas SW)")

fucus_BC <- plot_deriv(deriv1_fucus_BC, "deriv1 fucus BC")
fucus_NE <- plot_deriv(deriv1_fucus_NE, "deriv1 fucus NE")
fucus_NW <- plot_deriv(deriv1_fucus_NW, "deriv1 fucus NW")
fucus_SE <- plot_deriv(deriv1_fucus_SE, "deriv1 fucus SE")
fucus_SW <- plot_deriv(deriv1_fucus_SW, "deriv1 fucus SW")

littorina_BC <- plot_deriv(deriv1_littorina_BC, "deriv1 littorina BC")
littorina_NE <- plot_deriv(deriv1_littorina_NE, "deriv1 littorina NE")
littorina_NW <- plot_deriv(deriv1_littorina_NW, "deriv1 littorina NW")
littorina_SE <- plot_deriv(deriv1_littorina_SE, "deriv1 littorina SE")
littorina_SW <- plot_deriv(deriv1_littorina_SW, "deriv1 littorina SW")

mastocarpus_BC <- plot_deriv(deriv1_mastocarpus_BC, "deriv1 mastocarpus BC") ## Fix it
mastocarpus_NE <- plot_deriv(deriv1_mastocarpus_NE, "deriv1 mastocarpus NE")
mastocarpus_NW <- plot_deriv(deriv1_mastocarpus_NW, "deriv1 mastocarpus NW")
mastocarpus_SE <- plot_deriv(deriv1_mastocarpus_SE, "deriv1 mastocarpus SE")
mastocarpus_SW <- plot_deriv(deriv1_mastocarpus_SW, "deriv1 mastocarpus SW")

mytilus_BC <- plot_deriv(deriv1_mytilus_BC, "deriv1 mytilus BC")
#mytilus_MC <- plot_deriv(deriv1_mytilus_MC, "deriv1 mytilus MC")
mytilus_NE <- plot_deriv(deriv1_mytilus_NE, "deriv1 mytilus NE")
mytilus_NW <- plot_deriv(deriv1_mytilus_NW, "deriv1 mytilus NW")
mytilus_SE <- plot_deriv(deriv1_mytilus_SE, "deriv1 mytilus SE")
mytilus_SW <- plot_deriv(deriv1_mytilus_SW, "deriv1 mytilus SW")

nucella_BC <- plot_deriv(deriv1_nucella_BC, "deriv 1 nucella BC")
#nucella_MC <- plot_deriv(deriv1_nucella_MC, "deriv 2 nucella MC")
nucella_NE <- plot_deriv(deriv1_nucella_NE, "deriv 1 nucella NE")
nucella_NW <- plot_deriv(deriv1_nucella_NW, "deriv 1 nucella NW")
nucella_SE <- plot_deriv(deriv1_nucella_SE, "deriv 1 nucella SE")
nucella_SW <- plot_deriv(deriv1_nucella_SW, "deriv 1 nucella SW")

semibalanus_BC <- plot_deriv(deriv1_semibalanus_BC, "deriv 1 semibalanus BC")
semibalanus_NE <- plot_deriv(deriv1_semibalanus_NE, "deriv 1 semibalanus NE")
semibalanus_NW <- plot_deriv(deriv1_semibalanus_NW, "deriv 1 semibalanus NW")
semibalanus_SE <- plot_deriv(deriv1_semibalanus_SE, "deriv 1 semibalanus SE")
semibalanus_SW <- plot_deriv(deriv1_semibalanus_SW, "deriv 1 semibalanus SW")

amphipoda_BC <- plot_deriv(deriv1_amphipoda_BC, "deriv 1 amphipoda BC")
amphipoda_NE <- plot_deriv(deriv1_amphipoda_NE, "deriv 1 amphipoda NE")
amphipoda_NW <- plot_deriv(deriv1_amphipoda_NW, "deriv 1 amphipoda NW")
amphipoda_SE <- plot_deriv(deriv1_amphipoda_SE, "deriv 1 amphipoda SE")
amphipoda_SW <- plot_deriv(deriv1_amphipoda_SW, "deriv 1 amphipoda SW")

isopoda_BC <- plot_deriv(deriv1_isopoda_BC, "deriv 1 isopoda BC")
isopoda_NE <- plot_deriv(deriv1_isopoda_NE, "deriv 1 isopoda NE")
isopoda_NW <- plot_deriv(deriv1_isopoda_NW, "deriv 1 isopoda NW")
isopoda_SE <- plot_deriv(deriv1_isopoda_SE, "deriv 1 isopoda SE")
isopoda_SW <- plot_deriv(deriv1_isopoda_SW, "deriv 1 isopoda SW")


### Make multipanel figure for each to put in annexe

cover_deriv1_plot_global <- ggarrange(cover_MDS1_global, cover_MDS2_global, cover_MDS3_global,
                               nrow = 1, ncol = 3)
count_deriv1_plot_global <- ggarrange(count_MDS1_global, count_MDS2_global, count_MDS3_global,
                                      nrow = 1, ncol = 3)

count_deriv1_plot_site_MDS1 <- ggarrange(count_MDS1_BC,count_MDS1_NE,
                                   count_MDS1_NW, count_MDS1_SE,
                                   count_MDS1_SW,
                                   ncol = 3, nrow = 2)
count_deriv1_plot_site_MDS2 <- ggarrange(count_MDS2_BC,count_MDS2_NE,
                                         count_MDS2_NW, count_MDS2_SE,
                                         count_MDS2_SW,
                                         ncol = 3, nrow = 2)
count_deriv1_plot_site_MDS3 <- ggarrange(count_MDS3_BC,count_MDS3_NE,
                                         count_MDS3_NW, count_MDS3_SE,
                                         count_MDS3_SW,
                                         ncol = 3, nrow = 2)

cover_deriv1_plot_site_MDS1 <- ggarrange(cover_MDS1_BC,cover_MDS1_NE,
                                         cover_MDS1_NW, cover_MDS1_SE,
                                         cover_MDS1_SW,
                                         ncol = 3, nrow = 2)
cover_deriv1_plot_site_MDS2 <- ggarrange(cover_MDS2_BC,cover_MDS2_NE,
                                         cover_MDS2_NW,cover_MDS2_SE,
                                         cover_MDS2_SW,
                                         ncol = 3, nrow = 2)
cover_deriv1_plot_site_MDS3 <- ggarrange(cover_MDS3_BC,cover_MDS3_NE,
                                         cover_MDS3_NW,cover_MDS3_SE,
                                         cover_MDS3_SW,
                                         ncol = 3, nrow = 2)

deriv1_chondrus <- ggarrange(chondrus_BC,chondrus_NE,chondrus_NW,
                             chondrus_SE, chondrus_SW,
                             nrow = 2, ncol = 3)
deriv1_cmaenas <- ggarrange(cmaenas_BC,cmaenas_NE,cmaenas_NW,
                            cmaenas_SE, cmaenas_SW,
                            nrow = 2, ncol = 3)
deriv1_littorina <- ggarrange(littorina_BC,littorina_NE,littorina_NW,
                              littorina_SE,littorina_SW,
                              nrow = 2, ncol = 3)
deriv1_mastocarpus <- ggarrange(mastocarpus_BC,mastocarpus_NE,
                                mastocarpus_NW,mastocarpus_SE,mastocarpus_SW,
                                nrow = 2, ncol = 3)
deriv1_nucella <- ggarrange(nucella_BC,nucella_NE,nucella_NW,nucella_SE
                            ,nucella_SW,
                            nrow = 2, ncol = 3)
deriv1_mytilus <- ggarrange(mytilus_BC,mytilus_NE,mytilus_NW,
                            mytilus_SE,mytilus_SW,
                            nrow = 2, ncol = 3)
deriv1_ascophylum <- ggarrange(ascophyllum_BC,ascophyllum_NE,ascophyllum_NW,
                               ascophyllum_SE,ascophyllum_SW,
                               nrow = 2, ncol = 3)
deriv1_semibalanus <- ggarrange(semibalanus_BC,semibalanus_NE,semibalanus_NW,
                                semibalanus_SE,semibalanus_SW,
                                nrow = 2, ncol = 3)
deriv1_fucus <- ggarrange(fucus_BC,fucus_NE,fucus_NW,fucus_SE,fucus_SW,
                          nrow = 2, ncol = 3)

deriv1_amphipoda <- ggarrange(amphipoda_BC,amphipoda_NE,amphipoda_NW,amphipoda_SE,amphipoda_SW,
                              nrow = 2, ncol = 3)

deriv1_isopoda <- ggarrange(isopoda_BC,isopoda_NE,isopoda_NW,isopoda_SE,isopoda_SW,
                            nrow = 2, ncol = 3)


## Save fig of second deriv
path = "time_series_outputs/figures/supplement/deriv1/"

ggsave("deriv1_cover_global.tiff", path = path, 
       plot = cover_deriv1_plot_global, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

ggsave("deriv1_count_global.tiff", path = path, 
       plot = count_deriv1_plot_global, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

ggsave("deriv1_cover_site_MDS1.tiff", path = path, 
       plot = cover_deriv1_plot_site_MDS1, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_cover_site_MDS2.tiff", path = path, 
       plot = cover_deriv1_plot_site_MDS2, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_cover_site_MDS3.tiff", path = path, 
       plot = cover_deriv1_plot_site_MDS3, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_count_site_MDS1.tiff", path = path, 
       plot = count_deriv1_plot_site_MDS1, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_count_site_MDS2.tiff", path = path, 
       plot = count_deriv1_plot_site_MDS2, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_count_site_MDS3.tiff", path = path, 
       plot = count_deriv1_plot_site_MDS3, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

ggsave("deriv1_ascophyllum.tiff", path = path, 
       plot = deriv1_ascophylum, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_chondrus.tiff", path = path, 
       plot = deriv1_chondrus, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_cmaenas.tiff", path = path, 
       plot = deriv1_cmaenas, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_fucus.tiff", path = path, 
       plot = deriv1_fucus, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_littorina.tiff", path = path, 
       plot = deriv1_littorina, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_mastocarpus.tiff", path = path, 
       plot = deriv1_mastocarpus, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_mytilus.tiff", path = path, 
       plot = deriv1_mytilus, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_nucella.tiff", path = path, 
       plot = deriv1_nucella, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_semibalanus.tiff", path = path, 
       plot = deriv1_semibalanus, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

ggsave("deriv1_amphipoda.tiff", path = path, 
       plot = deriv1_amphipoda, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_isopoda.tiff", path = path, 
       plot = deriv1_isopoda, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")



### Do the same with the zi part of the model ##################################

ascophyllum_BC_zi <- plot_deriv(deriv1_ascophylum_BC_zi, "deriv1(ascophyllum BC zi)")
ascophyllum_NE_zi <- plot_deriv(deriv1_ascophylum_NE_zi, "deriv1(ascophyllum NE zi)")
ascophyllum_NW_zi <- plot_deriv(deriv1_ascophylum_BC_zi, "deriv1(ascophyllum NW zi)")
ascophyllum_SE_zi <- plot_deriv(deriv1_ascophylum_NE_zi, "deriv1(ascophyllum SE zi)")
ascophyllum_SW_zi <- plot_deriv(deriv1_ascophylum_BC_zi, "deriv1(ascophyllum SW zi)")

chondrus_BC_zi <- plot_deriv(deriv1_chondrus_BC_zi, "deriv1(chondrus BC zi)")
#chondrus_MC <- plot_deriv(deriv1_chondrus_MC, "deriv1(chondrus MC)")
chondrus_NE_zi <- plot_deriv(deriv1_chondrus_NE_zi, "deriv1(chondrus NE zi)")
chondrus_NW_zi <- plot_deriv(deriv1_chondrus_NW_zi, "deriv1(chondrus NW zi)")
chondrus_SE_zi <- plot_deriv(deriv1_chondrus_SE_zi, "deriv1(chondrus SE zi)")
chondrus_SW_zi <- plot_deriv(deriv1_chondrus_SW_zi, "deriv1(chondrus SW zi)")

cmaenas_BC_zi <- plot_deriv(deriv1_cmaenas_BC_zi, "deriv1(C maenas BC zi)")
#cmaenas_MC <- plot_deriv(deriv1_cmaenas_MC, "deriv1(C maenas MC)")
cmaenas_NE_zi <- plot_deriv(deriv1_cmaenas_NE_zi, "deriv1(C maenas NE zi)")
cmaenas_NW_zi <- plot_deriv(deriv1_cmaenas_NW_zi, "deriv1(C maenas NW zi)")
cmaenas_SE_zi <- plot_deriv(deriv1_cmaenas_SE_zi, "deriv1(C maenas SE zi)")
cmaenas_SW_zi <- plot_deriv(deriv1_cmaenas_SW_zi, "deriv1(C maenas SW zi)")

fucus_BC_zi <- plot_deriv(deriv1_fucus_BC_zi, "deriv1 fucus BC zi")
fucus_NE_zi <- plot_deriv(deriv1_fucus_NE_zi, "deriv1 fucus NE zi")
fucus_NW_zi <- plot_deriv(deriv1_fucus_NW_zi, "deriv1 fucus NW zi")
fucus_SE_zi <- plot_deriv(deriv1_fucus_SE_zi, "deriv1 fucus SE zi")
fucus_SW_zi <- plot_deriv(deriv1_fucus_SW_zi, "deriv1 fucus SW zi")

littorina_BC_zi <- plot_deriv(deriv1_littorina_BC_zi, "deriv1 littorina BC zi")
littorina_NE_zi <- plot_deriv(deriv1_littorina_NE_zi, "deriv1 littorina NE zi")
littorina_NW_zi <- plot_deriv(deriv1_littorina_NW_zi, "deriv1 littorina NW zi")
littorina_SE_zi <- plot_deriv(deriv1_littorina_SE_zi, "deriv1 littorina SE zi")
littorina_SW_zi <- plot_deriv(deriv1_littorina_SW_zi, "deriv1 littorina SW zi")

mastocarpus_BC_zi <- plot_deriv(deriv1_mastocarpus_BC_zi, "deriv1 mastocarpus BC zi") ## Fix it
mastocarpus_NE_zi <- plot_deriv(deriv1_mastocarpus_NE_zi, "deriv1 mastocarpus NE zi")
mastocarpus_NW_zi <- plot_deriv(deriv1_mastocarpus_NW_zi, "deriv1 mastocarpus NW zi")
mastocarpus_SE_zi <- plot_deriv(deriv1_mastocarpus_SE_zi, "deriv1 mastocarpus SE zi")
mastocarpus_SW_zi <- plot_deriv(deriv1_mastocarpus_SW_zi, "deriv1 mastocarpus SW zi")

mytilus_BC_zi <- plot_deriv(deriv1_mytilus_BC_zi, "deriv1 mytilus BC zi")
mytilus_NE_zi <- plot_deriv(deriv1_mytilus_NE_zi, "deriv1 mytilus NE zi")
mytilus_NW_zi <- plot_deriv(deriv1_mytilus_NW_zi, "deriv1 mytilus NW zi")
mytilus_SE_zi <- plot_deriv(deriv1_mytilus_SE_zi, "deriv1 mytilus SE zi")
mytilus_SW_zi <- plot_deriv(deriv1_mytilus_SW_zi, "deriv1 mytilus SW zi")

nucella_BC_zi <- plot_deriv(deriv1_nucella_BC_zi, "deriv 1 nucella BC zi")
nucella_NE_zi <- plot_deriv(deriv1_nucella_NE_zi, "deriv 1 nucella NE zi")
nucella_NW_zi <- plot_deriv(deriv1_nucella_NW_zi, "deriv 1 nucella NW zi")
nucella_SE_zi <- plot_deriv(deriv1_nucella_SE_zi, "deriv 1 nucella SE zi")
nucella_SW_zi <- plot_deriv(deriv1_nucella_SW_zi, "deriv 1 nucella SW zi")

semibalanus_BC_zi <- plot_deriv(deriv1_semibalanus_BC_zi, "deriv 1 semibalanus BC zi")
semibalanus_NE_zi <- plot_deriv(deriv1_semibalanus_NE_zi, "deriv 1 semibalanus NE zi")
semibalanus_NW_zi <- plot_deriv(deriv1_semibalanus_NW_zi, "deriv 1 semibalanus NW zi")
semibalanus_SE_zi <- plot_deriv(deriv1_semibalanus_SE_zi, "deriv 1 semibalanus SE zi")
semibalanus_SW_zi <- plot_deriv(deriv1_semibalanus_SW_zi, "deriv 1 semibalanus SW zi")

amphipoda_BC_zi <- plot_deriv(deriv1_amphipoda_BC_zi, "deriv 1 amphipoda BC zi")
amphipoda_NE_zi <- plot_deriv(deriv1_amphipoda_NE_zi, "deriv 1 amphipoda NE zi")
amphipoda_NW_zi <- plot_deriv(deriv1_amphipoda_NW_zi, "deriv 1 amphipoda NW zi")
amphipoda_SE_zi <- plot_deriv(deriv1_amphipoda_SE_zi, "deriv 1 amphipoda SE zi")
amphipoda_SW_zi <- plot_deriv(deriv1_amphipoda_SW_zi, "deriv 1 amphipoda SW zi")

isopoda_BC_zi <- plot_deriv(deriv1_isopoda_BC_zi, "deriv 1 isopoda BC zi")
isopoda_NE_zi <- plot_deriv(deriv1_isopoda_NE_zi, "deriv 1 isopoda NE zi")
isopoda_NW_zi <- plot_deriv(deriv1_isopoda_NW_zi, "deriv 1 isopoda NW zi")
isopoda_SE_zi <- plot_deriv(deriv1_isopoda_SE_zi, "deriv 1 isopoda SE zi")
isopoda_SW_zi <- plot_deriv(deriv1_isopoda_SW_zi, "deriv 1 isopoda SW zi")



### Make multipanel figure for each to put in annexe

deriv1_chondrus_zi <- ggarrange(chondrus_BC_zi,chondrus_NE_zi,chondrus_NW_zi,
                                chondrus_SE_zi, chondrus_SW_zi,
                                nrow = 2, ncol = 3)
deriv1_cmaenas_zi <- ggarrange(cmaenas_BC_zi,cmaenas_NE_zi,cmaenas_NW_zi,
                               cmaenas_SE_zi, cmaenas_SW_zi,
                               nrow = 2, ncol = 3)
deriv1_littorina_zi <- ggarrange(littorina_BC_zi,littorina_NE_zi,littorina_NW_zi,
                                 littorina_SE_zi,littorina_SW_zi,
                                 nrow = 2, ncol = 3)
deriv1_mastocarpus_zi <- ggarrange(mastocarpus_BC_zi,mastocarpus_NE_zi,
                                   mastocarpus_NW_zi,mastocarpus_SE_zi,mastocarpus_SW_zi,
                                   nrow = 2, ncol = 3)
deriv1_nucella_zi <- ggarrange(nucella_BC_zi,nucella_NE_zi,nucella_NW_zi,nucella_SE_zi
                               ,nucella_SW_zi,
                               nrow = 2, ncol = 3)
deriv1_mytilus_zi <- ggarrange(mytilus_BC_zi,mytilus_NE_zi,mytilus_NW_zi,
                               mytilus_SE_zi,mytilus_SW_zi,
                               nrow = 2, ncol = 3)
deriv1_ascophylum_zi <- ggarrange(ascophyllum_BC_zi,ascophyllum_NE_zi,ascophyllum_NW_zi,
                                  ascophyllum_SE_zi,ascophyllum_SW_zi,
                                  nrow = 2, ncol = 3)
deriv1_semibalanus_zi <- ggarrange(semibalanus_BC_zi,semibalanus_NE_zi,semibalanus_NW_zi,
                                   semibalanus_SE_zi,semibalanus_SW_zi,
                                   nrow = 2, ncol = 3)
deriv1_fucus_zi <- ggarrange(fucus_BC_zi,fucus_NE_zi,fucus_NW_zi,fucus_SE_zi,fucus_SW_zi,
                             nrow = 2, ncol = 3)

deriv1_amphipoda_zi <- ggarrange(amphipoda_BC_zi,amphipoda_NE_zi,amphipoda_NW_zi,amphipoda_SE_zi,
                                 amphipoda_SW_zi,
                                 nrow = 2, ncol = 3)

deriv1_isopoda_zi <- ggarrange(isopoda_BC_zi,isopoda_NE_zi,isopoda_NW_zi,isopoda_SE_zi,
                               isopoda_SW_zi,
                               nrow = 2, ncol = 3)


## Save fig of second deriv

ggsave("deriv1_ascophyllum_zi.tiff", path = path, 
       plot = deriv1_ascophylum_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_chondrus_zi.tiff", path = path, 
       plot = deriv1_chondrus_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_cmaenas_zi.tiff", path = path, 
       plot = deriv1_cmaenas_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_fucus_zi.tiff", path = path, 
       plot = deriv1_fucus_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_littorina_zi.tiff", path = path, 
       plot = deriv1_littorina_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_mastocarpus_zi.tiff", path = path, 
       plot = deriv1_mastocarpus_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_mytilus_zi.tiff", path = path, 
       plot = deriv1_mytilus_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_nucella_zi.tiff", path = path, 
       plot = deriv1_nucella_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_semibalanus_zi.tiff", path = path, 
       plot = deriv1_semibalanus_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

ggsave("deriv1_amphipoda_zi.tiff", path = path, 
       plot = deriv1_amphipoda_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
ggsave("deriv1_isopoda_zi.tiff", path = path, 
       plot = deriv1_isopoda_zi, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
