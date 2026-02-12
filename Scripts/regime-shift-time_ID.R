# This script is to calculate the year where each regime shift happenes (deriv max)
# with a 95% CI.
# Contributor: Julien Beaulieu

library(brms)
library(tidybayes)
library(ggpubr)
library(bayestestR)

# load models results
noaa <- readRDS("time_series_outputs/noaa_pcoa_GS_3.rds")

# load model results for the NOAA data but when selecting species 
# also detected in the keen data to test if the results are similar

noaa <- readRDS("time_series_outputs/noaa_pcoa-KEEN_GS_3.rds")

## Extract info from conditional_smooth========================================

### NOAA #######################################################################
#___________________________________________________________________________

msms_noaa <- conditional_smooths(noaa, spaghetti = T)

#extract data of the draws (spaghetti)
f <- attributes(msms_noaa$`mu_MDS1: s(year)`)$spaghetti
f <- f[,c("year", "estimate__", "sample__")]


f_w_noaa_MDS1 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
f_w_noaa_MDS1[is.nan(f_w_noaa_MDS1)] <- 0


f <- attributes(msms_noaa$`mu_MDS2: s(year)`)$spaghetti
f <- f[,c("year", "estimate__", "sample__")]
f_w_noaa_MDS2 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
f_w_noaa_MDS2[is.nan(f_w_noaa_MDS2)] <- 0

f <- attributes(msms_noaa$`mu_MDS3: s(year)`)$spaghetti
f <- f[,c("year", "estimate__", "sample__")]
f_w_noaa_MDS3 <- reshape2::dcast(f, year ~ sample__, value.var = "estimate__") %>%
  as_tibble()
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
f_w_noaa_MDS2[is.nan(f_w_noaa_MDS3)] <- 0

## create function to derivate of each draw and CI of the derivative ================================================

Deriv <- function(data){
  
  mty1 <- data[,-1]
  
  #calculate delta y 
  mty2 <- mty1[-1,]
  mty1 <- mty1[-nrow(data),]
  m_delta_y <- mty2 - mty1
  
  #calculate delta x
  x1 <- data$year[-nrow(data)]
  x2 <- data$year[-1]
  m_delta_x <- x2-x1
  
  #calculate derivative
  mt_deriv <- as.data.frame(matrix(ncol = ncol(mty1), nrow = nrow(m_delta_y)))
  for (i in c(1:ncol(m_delta_y))) {
    mt_deriv[,i] <- m_delta_y[,i]/m_delta_x
  }  
  mt_deriv$year <- (x1 + m_delta_x/2)
  return(mt_deriv)   
} 

### Calculate derivative

deriv_MDS1 <- Deriv(f_w_noaa_MDS1)
deriv_MDS2 <- Deriv(f_w_noaa_MDS2)
deriv_MDS3 <- Deriv(f_w_noaa_MDS3)


#divide MDS2 in two for each regime shift
deriv_MDS2_A <- filter(deriv_MDS3, year < 2000)
deriv_MDS2_B <- filter(deriv_MDS3, year > 2000)


# find the max with CI

find_time <- function(data){
  max <- NA
  years <- NA
  data2 <- data[, -which(names(data) %in% "year")]
  for (i in c(1:ncol(data2))) {
   max[i] <- max(data2[,i]) #change to max for 1990 and 2000 and min for 2010
  }
  for (i in c(1:ncol(data2))) {
    for (j in c(1:nrow(data2))) {
      if(data2[j,i] == max[i]){years[i] <- data$year[j]}

    }
    
  }
  return(years)
}



time_MDS1 <- find_time(deriv_MDS1)
max(time_MDS1)
min(time_MDS1)
mean(time_MDS1)
sd(time_MDS1)
quantile(time_MDS1, c(.05,.1,.25,.5,.75,.9,.95))

time_MDS2_A <- find_time(deriv_MDS2_A)
max(time_MDS2_A)
min(time_MDS2_A)
mean(time_MDS2_A)
sd(time_MDS2_A)
quantile(time_MDS2_A, c(.05,.1,.25,.5,.75,.9,.95))

time_MDS2_B <- find_time(deriv_MDS2_B)
max(time_MDS2_B)
min(time_MDS2_B)
mean(time_MDS2_B)
sd(time_MDS2_B)
quantile(time_MDS2_B, c(.05,.1,.25,.5,.75,.9,.95))


### plot the deriv #####################################################

# recalculate derivative to have a different output
source("Scripts/first_deriv_function.R")

deriv_MDS1 <- Deriv(f_w_noaa_MDS1)
deriv_MDS2 <- Deriv(f_w_noaa_MDS3)

# plot ###

source("Scripts/plot_deriv.R")

plot_MDS1 <- plot_deriv(deriv_MDS1, ylab = "Deriv 1 Subtidal PCoA 1", col = "darkgoldenrod")
plot_MDS2 <- plot_deriv(deriv_MDS2, ylab = "Deriv 1 Subtidal PCoA 2", col = "darkred")

plot_deriv_sub <- ggarrange(plot_MDS1,plot_MDS2,
                           nrow = 1, ncol = 2)

ggsave("time_series_outputs/figures/supplement/deriv1/noaa_pcoa_deriv1.tiff",
       plot = plot_deriv_sub, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
