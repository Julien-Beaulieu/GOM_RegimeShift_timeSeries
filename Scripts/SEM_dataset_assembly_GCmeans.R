#-----------------------------------------------------#
#-SEM dataset assembly--------------------------------#
#-Contributor: Nicole Knight--------------------------#
#-Date: April 8 2021----------------------------------#
#-----------------------------------------------------#

#### Guide ####

# 1. Load libraries
# 2. Assemble dataframe
# 3. Build test models 
# 4. Model checking and plotting

#Note that for each section in the guide, you can click on the small triangle left of the #### i. Header #### to collapse it

####-1. Load libraries-####

library(tidyverse)
library(janitor)
library(reshape2)
library(dagitty)
library(ggdag)
library(ggplot2)
library(brms)
library(bayesplot)
library(vegan)

####-2. Assemble dataframe-####


####-NOAA trawl data-####


#Notes on column headers:

# avg_individ_PUE = the average per-unit-effort number of individuals caught across all transects within 100 km of Appledore
# avg_tot_mass_PUE = the average per-unit-effort biomass of all individuals caught across all transects within 100 km of Appledore
# scaled_avg_individ_PUE = scaled (but not centered) avg_individ_PUE
# scaled_avg_tot_mass_PUE = scaled (but not centered) avg_tot_mass_PUE

combined_int_url <- "https://raw.githubusercontent.com/Intertidal-Subtidal-WG/additional_data_sources/master/data/per_season2.csv"
download.file(combined_int_url, destfile = "data/NOAA_trawl_data.csv")

noaa_data <- read_csv("data/NOAA_trawl_data.csv") %>%
  filter(scaled_med_individ_PUE != "NaN") %>% #there are some NaNs in the data but don't include any of our species of interest
  filter(scaled_med_tot_mass_PUE != "NaN") %>%
  mutate(scaled_med_individ_PUE = as.numeric(scaled_med_individ_PUE)) %>%
  mutate(scaled_med_tot_mass_PUE = as.numeric(scaled_med_tot_mass_PUE))

noaa_data_sub <- subset(noaa_data, SCINAME %in% c("CANCER BOREALIS", "CANCER IRRORATUS", "GADUS MORHUA", "TAUTOGOLABRUS ADSPERSUS", "HOMARUS AMERICANUS", "MYOXOCEPHALUS SCORPIUS", "POLLACHIUS VIRENS"))

noaa_data_sub <- aggregate(noaa_data_sub$scaled_med_tot_mass_PUE, by = list(noaa_data_sub$YEAR, noaa_data_sub$SCINAME), mean)

colnames(noaa_data_sub) <- c("YEAR", "SCINAME", "scaled_med_tot_mass_PUE")

#-Convert to wide dataframe-#
noaa_data_wide <- dcast(noaa_data_sub, YEAR ~ SCINAME, value.var = "scaled_med_tot_mass_PUE", mean) %>%
  as_tibble() %>%
  clean_names() %>%
  subset(year > 1981)

#-NaNs are years where no individuals of that species was detected, will replace with 0-#

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

noaa_data_wide[is.nan(noaa_data_wide)] <- 0


ggplot(noaa_data_wide, aes(x = year, y = gadus_morhua)) + geom_line(size = 2) + theme_bw(base_size = 16)
ggplot(noaa_data_wide, aes(x = year, y = cancer_borealis)) + geom_line(size = 2) + theme_bw(base_size = 16)
ggplot(noaa_data_wide, aes(x = year, y = tautogolabrus_adspersus)) + geom_line(size = 2) + theme_bw(base_size = 16)
ggplot(noaa_data_wide, aes(x = year, y = pollachius_virens)) + geom_line(size = 2) + theme_bw(base_size = 16)
ggplot(noaa_data_wide, aes(x = year, y = homarus_americanus)) + geom_line(size = 2) + theme_bw(base_size = 16)

noaa_data_wide <- noaa_data_wide %>%
  rename(cancer_irroratus_subtidal = cancer_irroratus)

noaa_data_wide$subtidal_crabs <- noaa_data_wide$cancer_borealis + noaa_data_wide$cancer_irroratus_subtidal

#-Download gull data-#

gull_dat <- read_csv("data/gull_data.csv") %>%
  clean_names()

colnames(gull_dat)[c(4:5)] <- c("total_gulls", "gull_notes")

gull_dat <- gull_dat[,1:4]

all_dat <- left_join(noaa_data_wide, gull_dat, by = "year")

all_dat$gulls_centered <- (all_dat$total_gulls - mean(all_dat$total_gulls, na.rm = TRUE))/sd(all_dat$total_gulls, na.rm = TRUE)


####-Temperature data-####

#-Download level-specific temperature data (note that right now levels don't match up with intertidal data)-#
#-Currently need to update the link each time, which is a real pain-#

#combined_int_url <- "https://raw.githubusercontent.com/Intertidal-Subtidal-WG/thermaltolerance/main/outputs/datasets/temps_shorelevels_bymonth_bysite.csv?token=ACIASFWUIQ6GYMYHHKC42G3A2T6B6"
#download.file(combined_int_url, destfile = "data/temps_shorelevels_bymonth_bysite.csv")

temp_dat <- read_csv("data/temps_shorelevels_bymonth_bysite.csv") %>%
  filter(metric == "mean") %>%
  filter(months > 5) %>%
  filter(months < 9)

temp_dat_agg <- aggregate(temp_dat$temp, by = list(temp_dat$years, temp_dat$shore_level, temp_dat$site), mean) %>%
  as_tibble()
colnames(temp_dat_agg) <- c("year", "level_m", "TEMP_site", "mean_summer_sst")

#-Create transect key to match back temperature site designations to intertidal site designations-#
#-No SE site in temperature data so will assign to NE Appledore for now-#

TEMP_site <- c("Northeast", "Northeast","Northeast", "Northwest", "Southwest", "Southwest")
site <- c("NE Appledore", "SE Appledore", "Malaga Cut", "NW Appledore", "SW Appledore", "Babb's Cove")

site_conv <- bind_cols(TEMP_site, site)
colnames(site_conv) <- c("TEMP_site", "site")

temp_dat_agg <- full_join(temp_dat_agg, site_conv, by = "TEMP_site")

#-Convert back level to the spatial scale of the original data-#

temp_dat_agg$level <- round((13 -(temp_dat_agg$level_m/0.348)), digits = 2)
(unique(temp_dat_agg$level))

ggplot(temp_dat_agg, aes(x = level, y = mean_summer_sst)) + geom_jitter(aes(colour = site)) + theme_bw() #looks right


#-Create group centered values for mean summer SST-#

sst_transect_means <- aggregate(temp_dat_agg$mean_summer_sst, by = list(temp_dat_agg$site), mean, na.rm = TRUE)
sst_transect_sds <- aggregate(temp_dat_agg$mean_summer_sst, by = list(temp_dat_agg$site), sd, na.rm = TRUE)
sst_transect <- full_join(sst_transect_means, sst_transect_sds, by = "Group.1")

colnames(sst_transect) <- c("site", "summer_sst_transect_mean", "summer_sst_transect_sd")

temp_dat_agg <- full_join(temp_dat_agg, sst_transect, by = "site")

temp_dat_agg$GC_sst <- (temp_dat_agg$mean_summer_sst - temp_dat_agg$summer_sst_transect_mean)/temp_dat_agg$summer_sst_transect_sd

all_dat <- full_join(all_dat, temp_dat_agg, by = c("year"))


#-Download aragonite saturation data from Salisbury & Jonsson 2018 (note that year_notlagged is the actual year the data were collected,-#
#-but we will be using aragonite saturation from the prior year as a predictor (e.g., using 1984 aragonite data to predict 1985 nucella data)-#
#-so the "year" column reflects that-#

arag_dat <- read_csv("data/aragonite_sat_salisbury2018.csv")
colnames(arag_dat)[1] <- "arag_data_collection_year"

all_dat <- full_join(all_dat, arag_dat, by = "year")

#-Other environmental data-#

env_data <- read_csv("data/env_data_annual_2sep20.csv") %>%
  clean_names()

env_data$year_collected <- env_data$year
env_data$year <- env_data$year_collected + 1 #will use temperature data from year before

all_dat <- full_join(all_dat, env_data, by = "year")


####-Subtidal habitat complexity data (from Djikstra et al. 2017)-####

cmplx_data <- read_csv("data/complexity_data.csv")

cmplx_data$YR <- cmplx_data$YEAR - min(cmplx_data$YEAR)

ggplot(cmplx_data, aes(x = YEAR, y = COMPLEXITY, colour = SITE)) + geom_point(aes(colour = SITE, shape = SITE), size = 6) + theme_bw(base_size = 16) + geom_smooth(method = "lm")

#Here I center/standardize the data to each site, then take the average for each year there is available data, then interpolate for the entire time period

app <- filter(cmplx_data, SITE == "Appledore Island")
app$COMPLEXITY_STD <- (app$COMPLEXITY - mean(app$COMPLEXITY))/sd(app$COMPLEXITY)

stari <- filter(cmplx_data, SITE == "Star Island")
stari$COMPLEXITY_STD <- (stari$COMPLEXITY - mean(stari$COMPLEXITY))/sd(stari$COMPLEXITY)

wht <- filter(cmplx_data, SITE == "White Island")
wht$COMPLEXITY_STD <- (wht$COMPLEXITY - mean(wht$COMPLEXITY))/sd(wht$COMPLEXITY)

cmplx_data <- bind_rows(app, stari, wht)

cmplx_ind <- aggregate(cmplx_data$COMPLEXITY_STD, by = list(cmplx_data$YEAR), mean)
colnames(cmplx_ind) <- c("YEAR", "COMPLEXITY_STD")

cmplx_int <- approx(cmplx_ind$YEAR, cmplx_ind$COMPLEXITY_STD, xout = 1977:2015, method = "linear")

cmplx_int <- as_tibble(cmplx_int)
colnames(cmplx_int) <- c("YEAR", "COMPLEXITY_INT")

cmplx_int <- left_join(cmplx_int, cmplx_ind, by = "YEAR")
colnames(cmplx_int) <- c("year", "complexity_int", "complexity_std")

ggplot(data = cmplx_int) + geom_point(aes(x = year, y = complexity_std), size = 5) + 
  geom_path(aes(x = year, y = complexity_int), size = 1, colour = "maroon") + theme_bw(base_size = 16) + ylab("Interpolated complexity index Djikstra et al. 2017")

all_dat <- full_join(all_dat, cmplx_int, by = "year")

####-Intertidal data-####


#-Note: I downloaded the file "combined_intertidal_abundance.RDS" directly from GitHub into my copy of the repo,-#
#-the script kept throwing errors if I used download.file()-#

intertidal_data_0 <- read_rds("data/combined_intertidal_abundance.RDS") 

intertidal_data <- filter(intertidal_data_0, ORGANISM %in% c("Semibalanus balanoides", "Mytilus edulis", "Nucella lapillus", 
                                                             "Littorina littorea", "Littorina obtusata", "Littorina saxatilis",
                                                             "Carcinus maenus", "Cancer irroratus", 
                                                             "Chondrus crispus",
                                                             "Ascophyllum nodosum", "Mastocarpus stellatus", 
                                                             "Fucus vesiculosus", "Fucus distichus"))  

#-Note that I'm not including Cancer borealis because it only appears 12 times in the entire dataset,-#
#-we will use NOAA trawl data for Cancer borealis-#

#-Only want percent cover data for mussels and barnacles-#

intertidal_data <- intertidal_data[!(intertidal_data$PROTOCOL == "Intertidal_Count" & intertidal_data$ORGANISM == "Semibalanus balanoides"), ]
intertidal_data <- intertidal_data[!(intertidal_data$PROTOCOL == "Intertidal_Count" & intertidal_data$ORGANISM == "Mytilus edulis"), ]

#-Get rid of impossible values-#

intertidal_data <- intertidal_data[!(intertidal_data$ORGANISM == "Semibalanus balanoides" & intertidal_data$VALUE > 100), ]
intertidal_data <- intertidal_data[!(intertidal_data$ORGANISM == "Mytilus edulis" & intertidal_data$VALUE > 100), ]


#-There are two plots with 30 and 60 C. maenas.  Could be an entry error or a cast of babies, but either way it's really messing with the estimates, so will remove-#

intertidal_data <- intertidal_data[!(intertidal_data$ORGANISM == "Carcinus maenus" & intertidal_data$VALUE > 29), ]


#-Recode 100% mussel or barnacle cover as 99% cover because beta distribution can't handle 1s (these data are converted to a 0-0.99 range later in the script)-#

length(intertidal_data$VALUE[which(intertidal_data$ORGANISM == "Mytilus edulis" & intertidal_data$VALUE == 100)]) #Recoding 11 entries
length(intertidal_data$VALUE[which(intertidal_data$ORGANISM == "Semibalanus balanoides" & intertidal_data$VALUE == 100)]) #Recoding 49 entries

intertidal_data$VALUE[which(intertidal_data$ORGANISM == "Mytilus edulis" & intertidal_data$VALUE == 100)] <- 99
intertidal_data$VALUE[which(intertidal_data$ORGANISM == "Semibalanus balanoides" & intertidal_data$VALUE == 100)] <- 99

#-Create a column with values group-centered (but not standardized) by transect and level-#

intertidal_means <- aggregate(intertidal_data$VALUE, by = list(intertidal_data$ORGANISM, intertidal_data$INTERTIDAL_TRANSECT, intertidal_data$LEVEL), mean)
colnames(intertidal_means) <- c("ORGANISM", "INTERTIDAL_TRANSECT", "LEVEL", 'GROUP_MEAN')

intertidal_data <- left_join(intertidal_data, intertidal_means, by = c("ORGANISM", "INTERTIDAL_TRANSECT", "LEVEL"))

#intertidal_data <- intertidal_data[!(intertidal_data$ORGANISM == "Semibalanus balanoides" & intertidal_data$VALUE > 100), ]
#intertidal_data <- intertidal_data[!(intertidal_data$ORGANISM == "Mytilus edulis" & intertidal_data$VALUE > 100), ]

intertidal_data_wide <- dcast(intertidal_data, SITE + INTERTIDAL_TRANSECT + YEAR + LEVEL + REPLICATE ~ ORGANISM, value.var = "VALUE", sum) %>%
  as_tibble() %>%
  clean_names() %>%
  rename(carcinus_maenas = carcinus_maenus) %>%
  rename(cancer_irroratus_intertidal = cancer_irroratus)

intertidal_data_wide_means <- dcast(intertidal_data, SITE + INTERTIDAL_TRANSECT + YEAR + LEVEL + REPLICATE ~ ORGANISM, value.var = "GROUP_MEAN", mean) %>%
  as_tibble() %>%
  clean_names() %>%
  rename(carcinus_maenas = carcinus_maenus) %>%
  rename(cancer_irroratus_intertidal = cancer_irroratus)

colnames(intertidal_data_wide_means)[6:18] <- paste(colnames(intertidal_data_wide_means)[6:18], "gc_mean", sep = "_")

#Convert percentages to fractions for zero-inflated beta distribution

intertidal_data_wide$mytilus_edulis <- intertidal_data_wide$mytilus_edulis/100
intertidal_data_wide$semibalanus_balanoides <- intertidal_data_wide$semibalanus_balanoides/100

intertidal_data_wide$ascophyllum_nodosum <- intertidal_data_wide$ascophyllum_nodosum/100
intertidal_data_wide$chondrus_crispus <- intertidal_data_wide$chondrus_crispus/100
intertidal_data_wide$fucus_distichus <- intertidal_data_wide$fucus_distichus/100
intertidal_data_wide$fucus_vesiculosus <- intertidal_data_wide$fucus_vesiculosus/100
intertidal_data_wide$mastocarpus_stellatus <- intertidal_data_wide$mastocarpus_stellatus/100


intertidal_data_wide_means$mytilus_edulis_gc_mean <- intertidal_data_wide_means$mytilus_edulis_gc_mean/100
intertidal_data_wide_means$semibalanus_balanoides_gc_mean <- intertidal_data_wide_means$semibalanus_balanoides_gc_mean/100

intertidal_data_wide_means$ascophyllum_nodosum_gc_mean <- intertidal_data_wide_means$ascophyllum_nodosum_gc_mean/100
intertidal_data_wide_means$chondrus_crispus_gc_mean <- intertidal_data_wide_means$chondrus_crispus_gc_mean/100
intertidal_data_wide_means$fucus_distichus_gc_mean <- intertidal_data_wide_means$fucus_distichus_gc_mean/100
intertidal_data_wide_means$fucus_vesiculosus_gc_mean <- intertidal_data_wide_means$fucus_vesiculosus_gc_mean/100
intertidal_data_wide_means$mastocarpus_stellatus_gc_mean <- intertidal_data_wide_means$mastocarpus_stellatus_gc_mean/100

intertidal_data_wide <- full_join(intertidal_data_wide, intertidal_data_wide_means, by = c("site", "intertidal_transect", "year", "level", "replicate"))

intertidal_data_wide$fucus <- intertidal_data_wide$fucus_distichus + intertidal_data_wide$fucus_vesiculosus
intertidal_data_wide$fucus_gc_mean <- intertidal_data_wide$fucus_distichus_gc_mean + intertidal_data_wide$fucus_vesiculosus_gc_mean

intertidal_data_wide$littorina <- intertidal_data_wide$littorina_littorea + intertidal_data_wide$littorina_obtusata + intertidal_data_wide$littorina_saxatilis
intertidal_data_wide$littorina_gc_mean <- intertidal_data_wide$littorina_littorea_gc_mean + intertidal_data_wide$littorina_obtusata_gc_mean + intertidal_data_wide$littorina_saxatilis_gc_mean

intertidal_data_wide$ascophyllum_nodosum_GC <- intertidal_data_wide$ascophyllum_nodosum - intertidal_data_wide$ascophyllum_nodosum_gc_mean
intertidal_data_wide$cancer_irroratus_GC <- intertidal_data_wide$cancer_irroratus - intertidal_data_wide$cancer_irroratus_gc_mean
intertidal_data_wide$carcinus_maenas_GC <- intertidal_data_wide$carcinus_maenas - intertidal_data_wide$carcinus_maenas_gc_mean
intertidal_data_wide$chondrus_crispus_GC <- intertidal_data_wide$chondrus_crispus - intertidal_data_wide$chondrus_crispus_gc_mean
intertidal_data_wide$fucus_distichus_GC <- intertidal_data_wide$fucus_distichus - intertidal_data_wide$fucus_distichus_gc_mean
intertidal_data_wide$fucus_vesiculos_GC <- intertidal_data_wide$fucus_vesiculos - intertidal_data_wide$fucus_vesiculos_gc_mean
intertidal_data_wide$littorina_littorea_GC <- intertidal_data_wide$littorina_littorea - intertidal_data_wide$littorina_littorea_gc_mean
intertidal_data_wide$littorina_obtusata_GC <- intertidal_data_wide$littorina_obtusata - intertidal_data_wide$littorina_obtusata_gc_mean
intertidal_data_wide$littorina_saxatilis_GC <- intertidal_data_wide$littorina_saxatilis - intertidal_data_wide$littorina_saxatilis_gc_mean
intertidal_data_wide$mastocarpus_stellatus_GC <- intertidal_data_wide$mastocarpus_stellatus - intertidal_data_wide$mastocarpus_stellatus_gc_mean
intertidal_data_wide$mytilus_edulis_GC <- intertidal_data_wide$mytilus_edulis - intertidal_data_wide$mytilus_edulis_gc_mean
intertidal_data_wide$nucella_lapillus_GC <- intertidal_data_wide$nucella_lapillus - intertidal_data_wide$nucella_lapillus_gc_mean
intertidal_data_wide$semibalanus_balanoides_GC <- intertidal_data_wide$semibalanus_balanoides - intertidal_data_wide$semibalanus_balanoides_gc_mean
intertidal_data_wide$fucus_GC <- intertidal_data_wide$fucus - intertidal_data_wide$fucus_gc_mean
intertidal_data_wide$littorina_GC <- intertidal_data_wide$littorina - intertidal_data_wide$littorina_gc_mean

#-Join intertidal data with other datasets-#

all_dat <- full_join(intertidal_data_wide, all_dat, by = c("year", "site", "level"))


###-More data transformations-####

#-Divide level predictor by 10 to make model fitting easier-#

all_dat$level_predictor <- all_dat$level/10

#-Create an "exposed" variable to distinguish sites on East/West Appledore-#

exposed <- c(1, 1, 1, 0, 0, 0)
site <- c("NE Appledore", "SE Appledore", "Malaga Cut", "NW Appledore", "SW Appledore", "Babb's Cove")

exp_table <- bind_cols(site, exposed)
colnames(exp_table) <- c("site", "exposed")

all_dat <- full_join(all_dat, exp_table, by = "site")

#-Center and standardize summer SST data-#

all_dat$sst_summer_std <- (all_dat$sst_summer - mean(all_dat$sst_summer, na.rm = TRUE))/sd(all_dat$sst_summer, na.rm = TRUE)

#-Clean out some impossible values-#

all_dat$mollusc_cover <- all_dat$mytilus_edulis + all_dat$semibalanus_balanoides #Clean out some impossible values

all_dat <- filter(all_dat, mollusc_cover < 1.001)

write.csv(all_dat, file = "data/full_dataset_for_SEMs_groupcentered.csv")
