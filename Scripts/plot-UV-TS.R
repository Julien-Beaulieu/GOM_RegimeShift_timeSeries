## This script is to code graph for species-specific time series
## The synchrony was determined by the script second_derivative_central_diff.R
# contributor: Julien Beaulieu

#load packages

library(brms)
library(ggplot2)
library(grid)
library(ggpubr)
library(janitor)
library(data.table)
library(boot)
library(dplyr)
library(tidyr)
library(stringr)

#load models outputs
# Best models

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


# load table with deriv 2 values

deriv_values <- read.csv("time_series_outputs/values_deriv2_2.csv")
deriv <- deriv_values %>% select(-X) %>% 
  mutate(specie = 
           ifelse(str_detect(var,"Ascophyllum"),"Ascophyllum",
                  ifelse(str_detect(var,"maenas"),"C. maenas",
                         ifelse(str_detect(var,"Chondrus"),"Chondrus",
                                ifelse(str_detect(var,"Fucus"),"Fucus",
                                       ifelse(str_detect(var,"Littorina"),"Littorina",
                                              ifelse(str_detect(var,"Mastocarpus"),"Mastocarpus",
                                                     ifelse(str_detect(var,"ytilus"),"Mytilus",
                                                            ifelse(str_detect(var,"Nucella"),"Nucella",
                                                                          ifelse(str_detect(var,"Semibalanus"),"Semibalanus",
                                                                                 ifelse(str_detect(var,"Amphipoda"),"Amphipoda",
                                                                                        ifelse(str_detect(var,"Isopoda"),"Isopoda","-")))))))))))) %>% 
  mutate(site = ifelse(str_detect(var,"NW"),"NW Appledore",
                       ifelse(str_detect(var,"SW"),"SW Appledore",
                              ifelse(str_detect(var,"NE"),"NE Appledore",
                                     ifelse(str_detect(var,"SE"),"SE Appledore","Babb's Cove")))))

# calculate marginal pp
msms_ascophyllum <- conditional_smooths(ascophyllum)
msms_chondrus <- conditional_smooths(chondrus)
msms_cmenas <- conditional_smooths(cmenas)
msms_fucus <- conditional_smooths(fucus)
msms_littorina <- conditional_smooths(littorina)
msms_mastocarpus <- conditional_smooths(mastocarpus)
msms_mytilus <- conditional_smooths(mytilus)
msms_nucella <- conditional_smooths(nucella)
msms_semibalanus <- conditional_smooths(semibalanus)
msms_amphipoda <- conditional_smooths(amphipoda)
msms_isopoda <- conditional_smooths(isopoda)

# extract pp data for each model

ascophylum_pp <- msms_ascophyllum$`mu: s(year,by=site,bs="tp",m=2)`
ascophylum_pp_zi <- msms_ascophyllum$`hu: s(year,by=site,bs="tp",m=2)`

chondrus_pp <- msms_chondrus$`mu: s(year,by=site,bs="tp",m=2)`
chondrus_pp_zi <- msms_chondrus$`hu: s(year,by=site,bs="tp",m=2)`

cmenas_pp <- msms_cmenas$`mu: s(year,by=site,bs="tp",m=2)`
cmenas_pp_zi <- msms_cmenas$`zi: s(year,by=site,bs="tp",m=2)`

fucus_pp <- msms_fucus$`mu: s(year,by=site,bs="tp",m=2)`
fucus_pp_zi <- msms_fucus$`hu: s(year,by=site,bs="tp",m=2)`

littorina_pp <- msms_littorina$`mu: s(year,by=site,bs="tp",m=2)`
littorina_pp_zi <- msms_littorina$`zi: s(year,by=site,bs="tp",m=2)`

mastocarpus_pp <- msms_mastocarpus$`mu: s(year,by=site,bs="tp",m=2)`
mastocarpus_pp_zi <- msms_mastocarpus$`hu: s(year,by=site,bs="tp",m=2)`

mytilus_pp <- msms_mytilus$`mu: s(year,by=site,bs="tp",m=2)`
mytilus_pp_zi <- msms_mytilus$`zi: s(year,by=site,bs="tp",m=2)`

nucella_pp <- msms_nucella$`mu: s(year,by=site,bs="tp",m=2)`
nucella_pp_zi <- msms_nucella$`zi: s(year,by=site,bs="tp",m=2)`

semibalanus_pp <- msms_semibalanus$`mu: s(year,by=site,bs="tp",m=2)`
semibalanus_pp_zi <- msms_semibalanus$`zi: s(year,by=site,bs="tp",m=2)`

amphipoda_pp <- msms_amphipoda$`mu: s(year,by=site,bs="tp",m=2)`
amphipoda_pp_zi <- msms_amphipoda$`zi: s(year,by=site,bs="tp",m=2)`

isopoda_pp <- msms_isopoda$`mu: s(year,by=site,bs="tp",m=2)`
isopoda_pp_zi <- msms_isopoda$`zi: s(year,by=site,bs="tp",m=2)`


### remove MC before re-running
amphipoda_pp <- amphipoda_pp %>% filter(site != "Malaga Cut")
amphipoda_pp_zi <- amphipoda_pp_zi %>% filter(site != "Malaga Cut")

isopoda_pp <- isopoda_pp %>% filter(site != "Malaga Cut")
isopoda_pp_zi <- isopoda_pp_zi %>% filter(site != "Malaga Cut")

## Back transform y axis ######################################################

back_transform <- function(data, link_function){
  if(link_function == "logit"){
      data$estimate__ <- inv.logit(data$estimate__)
      data$lower__ <- inv.logit(data$lower__)
      data$upper__ <- inv.logit(data$upper__)
      data$se__ <- inv.logit(data$se__)
  }
  if(link_function == "log"){
    data$estimate__ <- exp(data$estimate__)
    data$lower__ <- exp(data$lower__)
    data$upper__ <- exp(data$upper__)
    data$se__ <- exp(data$se__)
  }
  return(data)
}
# 
# ascophylum_pp <- back_transform(data = ascophylum_pp, link_function = "log")
# chondrus_pp <- back_transform(chondrus_pp, link_function = "log")
# chondrus_pp_site <- back_transform(chondrus_pp_site, link_function = "log")
# cmenas_pp <- back_transform(cmenas_pp, link_function = "log")
# fucus_pp <- back_transform(fucus_pp, link_function = "log")
# littorina_pp <- back_transform(littorina_pp,"log")
# mastocarpus_pp <- back_transform(mastocarpus_pp, "log")
# mastocarpus_pp_site <- back_transform(mastocarpus_pp_site, "log")
# mytilus_pp <- back_transform(mytilus_pp, "logit")
# nucella_pp <- back_transform(nucella_pp, "log")
# semibalanus_pp <- back_transform(semibalanus_pp, "logit")


### Make function to plot #################################

plot_UV_TS <- function(pp = ascophylum_pp, der = deriv,sp = "Ascophyllum", 
                       type = "mu", save = TRUE){
  
  list_site <- unique(pp$site)
  
  for (i in 1:length(list_site)) {
    
    pp_site <- pp %>% filter(site == list_site[i])
    
    # make a color vector with grey if the 95% CI 2nd derivative overlap 0 
    # during the regime shift, and red if not
    col <- c("darkgrey","darkgrey","darkgrey")
    der_site <- der %>% filter(site == pp_site$site[1] & specie == sp)
    if(type == "mu"){
      for (j in c(1:nrow(der_site))) {
        if(der_site$upper[j] * der_site$lower[j] > 0){col[j] = "red"}
      }
    }
    if(type == "zi"){
      for (j in c(1:nrow(der_site))) {
        if(der_site$upper_zi[j] * der_site$lower_zi[j] > 0){col[j] = "red"}
      }
    }
    
    plot <- ggplot() +
      geom_line(data = pp_site,
                aes(y = estimate__, x = year), color = "darkblue",
                linewidth = 1) +
      geom_ribbon(data = pp_site, aes(
          y = estimate__,
          x = year,
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
      ylab(paste(type,sp, sep = " ")) +
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
    assign(paste(list_site[i]),plot)
  }
  
  mult_plot <- ggarrange(`Babb's Cove`,`NE Appledore`,`NW Appledore`,
                         `SE Appledore`, `SW Appledore`,
                         labels = c("Babb's Cove","NE Appledore",
                                    "NW Appledore","SE Appledore",
                                    "SW Appledore"),
                         hjust = -0.2, vjust = -0.075,
                         font.label = list(size = 15, face = "bold.italic"),
                         nrow = 2, ncol = 3)
  
  fig <- annotate_figure(mult_plot, top = text_grob("", color = "white", size = 30), 
  )
  fig
  
  if(save == TRUE){
    ggsave(paste("time_series_outputs/figures/supplement/UV_ts/UV_TS_",type,sp,".tiff", sep = ""),
           plot = fig, device = "tiff", scale = 2, 
           width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
           compression = "lzw")
  }
  
  return(fig)
}


### run the function to plot and save the figures #########################


plot_ascophyllum <- plot_UV_TS(save = TRUE)
plot_ascophyllum_zi <- plot_UV_TS(pp = ascophylum_pp_zi,
                                  type = "zi", save = TRUE)

plot_chondrus <- plot_UV_TS(pp = chondrus_pp, sp = "Chondrus", type = "mu")
plot_chondrus_zi <- plot_UV_TS(pp = chondrus_pp_zi, sp = "Chondrus", type = "zi")

plot_cmaenas <- plot_UV_TS(pp = cmenas_pp, sp = "C. maenas", type = "mu")
plot_cmaenas_zi <- plot_UV_TS(pp = cmenas_pp_zi, sp = "C. maenas", type = "zi")

plot_fucus <- plot_UV_TS(pp = fucus_pp, sp = "Fucus", type = "mu")
plot_fucus_zi <- plot_UV_TS(pp = fucus_pp_zi, sp = "Fucus", type = "zi")

plot_littorina <- plot_UV_TS(pp = littorina_pp, sp = "Littorina", type = "mu")
plot_littorina_zi <- plot_UV_TS(pp = littorina_pp_zi, sp = "Littorina", type = "zi")

plot_mastocarpus <- plot_UV_TS(pp = mastocarpus_pp, sp = "Mastocarpus", type = "mu")
plot_mastocarpus_zi <- plot_UV_TS(pp = mastocarpus_pp_zi, sp = "Mastocarpus", type = "zi")

plot_mytilus <- plot_UV_TS(pp = mytilus_pp, sp = "Mytilus", type = "mu")
plot_mytilus_zi <- plot_UV_TS(pp = mytilus_pp_zi, sp = "Mytilus", type = "zi")

plot_nucella <- plot_UV_TS(pp = nucella_pp, sp = "Nucella", type = "mu", save = FALSE)
plot_nucella_zi <- plot_UV_TS(pp = nucella_pp_zi, sp = "Nucella", type = "zi")

plot_semibalanus <- plot_UV_TS(pp = semibalanus_pp, sp = "Semibalanus", type = "mu")
plot_semibalanus_zi <- plot_UV_TS(pp = semibalanus_pp_zi, sp = "Semibalanus", type = "zi")

plot_amphipoda <- plot_UV_TS(pp = amphipoda_pp, sp = "Amphipoda", type = "mu")
plot_amphipoda_zi <- plot_UV_TS(pp = amphipoda_pp_zi, sp = "Amphipoda", type = "zi")

plot_isopoda <- plot_UV_TS(pp = isopoda_pp, sp = "Isopoda", type = "mu")
plot_isopoda_zi <- plot_UV_TS(pp = isopoda_pp_zi, sp = "Isopoda", type = "zi")

