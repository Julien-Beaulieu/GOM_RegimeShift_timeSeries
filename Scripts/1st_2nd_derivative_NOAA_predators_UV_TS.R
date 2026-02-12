# This script is to calculate the first derivative of the time series
# of subtidal predators directly feeding on intertidal preys. These 
# species were identified in dataPrep_NOAA_intertidal_Web.R and 
# NOAA-intertidal_interaction_web.R. The timeseries were fitted in 
# fit_NOAA_UV_TS_linkers.R.
# contributor: Julien Beaulieu

### load packages
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
library(bayestestR)
library(ggrepel)
source("Scripts/format_sp_name_for_globi.R")
source("Scripts/first_deriv_function.R")
source("Scripts/function_get_deriv_value.R")
source("Scripts/second_deriv_function.R")

### load data

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

###__________________________________________________________________________
### make a loop to load model output, plot the TS, calculate 1st derivative # 
### deriv2, and plot it and save the three figures and the driv values #######
###______________________________________________________________________

# make df for values

all_val_der_1 <- as.data.frame(matrix(nrow = 1,ncol = 4)) %>% 
  rename("var" = V1, "regime_shift" = V2, "upper" = V3, "lower" = V4)

all_val_der_2 <- as.data.frame(matrix(nrow = 1,ncol = 4)) %>% 
  rename("var" = V1, "regime_shift" = V2, "upper" = V3, "lower" = V4)


for (i in c(1:length(sp_list))) {
  # load model
  mod <- readRDS(paste("time_series_outputs/UV_TS_subtidal/",sp_list[i],".RDS", sep = ""))
  
  
  # get posterior predictions
  msms <- conditional_smooths(mod, spaghetti = TRUE, ndraws = 1000)
  PP <- msms$`mu: s(YEAR)`
  
  
  # calculate first derivative
  
  deriv1 <- attributes(msms$`mu: s(YEAR)`)$spaghetti %>% 
    dplyr::select(-c(effect1__, cond__)) %>% 
    pivot_wider(names_from = sample__, values_from = estimate__) %>% 
    rename("year" = YEAR) %>% 
    Deriv()
    
  # determine in which RS the 1st deriv is different from zero
  
  val_der1 <- values_deriv2(deriv1, name = sp_list[i])
  all_val_der_1 <- rbind(all_val_der_1, val_der1)
  col <- c("darkgrey","darkgrey","darkgrey")
  for (j in c(1:3)) {
    if(val_der1$upper[j] * val_der1$lower[j] > 0){col[j] = "red"}
  }
  
  # plot and save first derivative
  plot_der1 <- ggplot() +
    geom_line(data = deriv1,
              aes(y = y_deriv_estim, x = year), color = "darkblue",
              linewidth = 1) +
    geom_ribbon(data = deriv1, aes(
      y = y_deriv_estim,
      x = year,
      ymin = lower_CI,
      ymax = upper_CI
    ),
    fill = "darkblue",
    color = "darkblue",
    alpha = 0.15
    ) +
    geom_rect(aes(
      xmin = c(1987, 1999.4, 2009.5),
      xmax = c(1992.7, 2002, 2011.7),
      ymin = -Inf,
      ymax = Inf
    ),
    fill = col,
    alpha = 0.3) +
    ylab(paste("Deriv1", sp_list[i], sep = " ")) +
    xlab("Time") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      #panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),
      #remove minor-grid labels
      plot.background = element_blank(),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13)
    )
  
  ggsave(paste("time_series_outputs/UV_TS_subtidal/deriv1_trend_sub_UV_T/",
               sp_list[i],".tiff", sep=""),
         plot = plot_der1, device = "tiff", scale = 2, 
         width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
         compression = "lzw")
  
  
  #### calculate deriv 2 ########################
  
  deriv2 <- attributes(msms$`mu: s(YEAR)`)$spaghetti %>% 
    dplyr::select(-c(effect1__, cond__)) %>% 
    pivot_wider(names_from = sample__, values_from = estimate__) %>% 
    rename("year" = YEAR) %>% 
    Deriv2()
  
  # determine in which RS the 1st deriv is different from zero
  
  val_der2 <- values_deriv2(deriv2, name = sp_list[i])
  all_val_der_2 <- rbind(all_val_der_2,val_der2)
  col <- c("darkgrey","darkgrey","darkgrey")
  for (j in c(1:3)) {
    if(val_der2$upper[j] * val_der2$lower[j] > 0){col[j] = "red"}
  }
  
  
  # plot and save second derivative
  plot_der2 <- ggplot() +
    geom_line(data = deriv2,
              aes(y = y_deriv_estim, x = year), color = "darkblue",
              linewidth = 1) +
    geom_ribbon(data = deriv2, aes(
      y = y_deriv_estim,
      x = year,
      ymin = lower_CI,
      ymax = upper_CI
    ),
    fill = "darkblue",
    color = "darkblue",
    alpha = 0.15
    ) +
    geom_rect(aes(
      xmin = c(1987, 1999.4, 2009.5),
      xmax = c(1992.7, 2002, 2011.7),
      ymin = -Inf,
      ymax = Inf
    ),
    fill = col,
    alpha = 0.3) +
    ylab(paste("Deriv1", sp_list[i], sep = " ")) +
    xlab("Time") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      #panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),
      #remove minor-grid labels
      plot.background = element_blank(),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13)
    )
  
  ggsave(paste("time_series_outputs/UV_TS_subtidal/der2_trend_sub_UV_TS/",
               sp_list[i],".tiff", sep=""),
         plot = plot_der2, device = "tiff", scale = 2, 
         width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
         compression = "lzw")
  
  
  
  #### plot and save temporal trend ######
  
  
  plot <- ggplot() +
    geom_line(data = PP,
              aes(y = estimate__, x = YEAR), color = "darkblue",
              linewidth = 1) +
    geom_ribbon(data = PP, aes(
      y = estimate__,
      x = YEAR,
      ymin = lower__,
      ymax = upper__
    ),
    fill = "darkblue",
    color = "darkblue",
    alpha = 0.15
    ) +
    geom_rect(aes(
      xmin = c(1987, 1999.4, 2009.5),
      xmax = c(1992.7, 2002, 2011.7),
      ymin = -Inf,
      ymax = Inf
    ),
    fill = col,
    alpha = 0.3) +
    ylab(paste("mu", sp_list[i], sep = " ")) +
    xlab("Time") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      #panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),
      #remove minor-grid labels
      plot.background = element_blank(),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13)
    )
  
  ggsave(paste("time_series_outputs/UV_TS_subtidal/TS_trend_sub_UV_TS/",
               sp_list[i],".tiff", sep=""),
         plot = plot, device = "tiff", scale = 2, 
         width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
         compression = "lzw")
  print(i)
  print(sp_list[i])
  
  
}

all_val_der_1 <- all_val_der_1[-1,]

all_val_der_2 <- all_val_der_2[-1,]

write.csv(all_val_der_1, 
          "time_series_outputs/UV_TS_subtidal/deriv1_trend_sub_UV_T/der1_values_sub_UV_TS.csv")
write.csv(all_val_der_2,
          "time_series_outputs/UV_TS_subtidal/der2_trend_sub_UV_TS/der2_values_sub_UV_TS.csv")
