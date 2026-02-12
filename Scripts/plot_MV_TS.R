# This script si to plot time series on ordination axis.
# The models used are those selected in model_selection_mvTS.R
# the 2nd derivative values are from
# contributor: Julien Beaulieu

# Load packages

library(ggpubr)
library(grid)
library(ggplot2)
library(tidyverse)
library(rstan)
library(devtools)
library(brms)
library(rstan)
library(rstantools)
library(schoenberg)
library(bayesplot)
library(data.table)
library(reshape2)
library(bayestestR)
library(geomtextpath)
library(tidybayes)
library(wesanderson)

# Load model outputs

# load models results
# cover_GS <- readRDS("time_series_outputs/new/cover_pcoa_GS_9.rds")
# count_4 <- readRDS("time_series_outputs/new/count_pcoa_mix_9.rds")
count <- readRDS("time_series_outputs/output/count_pcoa_GI_10.rds")
cover <- readRDS("time_series_outputs/output/cover_pcoa_GI_10.rds")
noaa <- readRDS("time_series_outputs/output/noaa_pcoa_GS_3.rds")

# noaa when taking only the species present in the Keen data
noaa <- readRDS("time_series_outputs/noaa_pcoa-KEEN_GS_3.rds")

# load deriv 2 values for colors

deriv_values <- read.csv("time_series_outputs/values_deriv2_2.csv")
deriv <- deriv_values %>% select(-X) %>% 
  mutate(specie = 
           ifelse(str_detect(var,"Sessile MDS1"),"Sessile MDS1",
                  ifelse(str_detect(var,"Sessile MDS2"),"Sessile MDS2",
                         ifelse(str_detect(var,"Sessile MDS3"),"Sessile MDS3",
                                ifelse(str_detect(var,"Mobile MDS1"),"Mobile MDS1",
                                       ifelse(str_detect(var,"Mobile MDS2"),"Mobile MDS2",
                                              ifelse(str_detect(var,"Mobile MDS3"),"Mobile MDS3","-"))))))) %>% 
  mutate(site = ifelse(str_detect(var,"NW"),"NW Appledore",
                       ifelse(str_detect(var,"SW"),"SW Appledore",
                              ifelse(str_detect(var,"NE"),"NE Appledore",
                                     ifelse(str_detect(var,"SE"),"SE Appledore",
                                            ifelse(str_detect(var,"BC"), "Babb's Cove", "global"))))))

## Extract posterior prediction ===========================================

#extract data of the draws (spaghetti)

msms_cover <- conditional_smooths(cover)
msms_count <- conditional_smooths(count)

### cover global trend ###############################################

cover_MDS1_global <- msms_cover$`mu_MDS1: s(YEAR,bs="tp")`
cover_MDS2_global <- msms_cover$`mu_MDS2: s(YEAR,bs="tp")`
cover_MDS3_global <- msms_cover$`mu_MDS3: s(YEAR,bs="tp`

### count global trend ############################################

count_MDS1_global <- msms_count$`mu_MDS1: s(YEAR,bs="tp",m=3)`
count_MDS2_global <- msms_count$`mu_MDS2: s(YEAR,bs="tp",m=3)`
count_MDS3_global <- msms_count$`mu_MDS3: s(YEAR,bs="tp",m=3)`


### cover site specific #####################################

### Function to extract pp for univariate models

site_cover_MDS1 <- msms_cover$`mu_MDS1: s(YEAR,by=SITE,m=1,bs="tp")`
site_cover_MDS2 <- msms_cover$`mu_MDS2: s(YEAR,by=SITE,m=1,bs="tp")`
site_cover_MDS3 <- msms_cover$`mu_MDS3: s(YEAR,by=SITE,m=1,bs="tp")`


### Count site specific ####################################################

site_count_MDS1 <- msms_count$`mu_MDS1: s(YEAR,by=SITE,m=1,bs="tp")`
site_count_MDS2 <- msms_count$`mu_MDS2: s(YEAR,by=SITE,m=1,bs="tp")`
site_count_MDS3 <- msms_count$`mu_MDS3: s(YEAR,by=SITE,m=1,bs="tp")`

## NOAA = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

msms_noaa <- conditional_smooths(noaa, spaghetti = T)
noaa_MDS1 <- msms_noaa$`mu_MDS1: s(year)`
noaa_MDS2 <- msms_noaa$`mu_MDS2: s(year)`

### Plot ======================================================================

# Create dummy data to plot arrows with taxa corelated with axis

dum_data <- as.data.frame(cbind(x=c(1,1,2,2), y=c(1,8.5,1,8.5), axis=c("intertidal","intertidal","subtidal","subtidal")))
dum_data$x <- as.numeric(dum_data$x)
dum_data$y <- as.numeric(dum_data$y)



# make function to plot 2nd derivative for the intertidal



plot_2 <- function(pp, sp, der = deriv,
                   y_lab = paste("mu", sp, sep = " "),
                   save = TRUE, filename) {
  
  # Color by group
  main_col <- ifelse(str_detect(sp, "Mobile"), "#00A08A", "#046C9A")
  
  # Check if SITE exists
  if(!("SITE" %in% colnames(pp))) {
    pp$SITE <- "global"
  }
  
  pp <- pp %>% rename("site" = SITE)
  list_site <- unique(pp$site)
  
  plot_list <- list()
  
  for (i in seq_along(list_site)) {
    pp_site <- pp %>% filter(site == list_site[i])
    
    # Color for rectangles
    colvect <- rep("darkgrey", 3)
    der_site <- der %>% filter(site == pp_site$site[1], specie == sp) %>% unique()
    for (j in 1:nrow(der_site)) {
      if (der_site$upper[j] * der_site$lower[j] > 0) {
        colvect[j] <- "#FF0000"
      }
    }
    
    rect_df <- data.frame(
      xmin = c(1987, 1999.4, 2009.5),
      xmax = c(1992.7, 2002, 2011.7),
      ymin = -Inf,
      ymax = Inf,
      colvect = colvect
    )
    
    p <- ggplot() +
      geom_line(data = pp_site,
                aes(y = estimate__, x = YEAR),
                color = main_col, linewidth = 1) +
      geom_ribbon(data = pp_site,
                  aes(y = estimate__, x = YEAR, ymin = lower__, ymax = upper__),
                  fill = main_col, color = main_col, alpha = 0.25) +
      geom_rect(data = rect_df,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colvect),
                alpha = 0.25, inherit.aes = FALSE) +
      scale_fill_identity() +
      ylab(y_lab) +
      xlab("Time") +
      geom_hline(yintercept = 0) +
      theme_bw() +
      theme(
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.y = element_blank()
      )
    
    plot_list[[list_site[i]]] <- p
  }
  
  # Combine plots
  if ("global" %in% list_site) {
    final_plot <- plot_list[["global"]]  # just the ggplot object
  } else {
    final_plot <- ggarrange(plotlist = plot_list,
                            labels = list_site,
                            font.label = list(size = 15, face = "bold.italic"),
                            nrow = 2, ncol = 3)
  }
  
  # Save if requested
  if(save && !is.null(filename)) {
    ggsave(paste0("time_series_outputs/figures/", filename, ".tiff"),
           plot = final_plot, device = "tiff", scale = 2,
           width = 1920, height = 991, units = "px", dpi = 300,
           bg = "white", compression = "lzw")
  }
  
  return(final_plot)
}



#### use the function to plot and save figures  

plot_count_MDS1_glo <- plot_2(pp = count_MDS1_global, sp = "Mobile MDS1", save = FALSE,
                              y_lab = "mu Mobile PCoA1 (52%)")
plot_count_MDS2_glo <- plot_2(pp = count_MDS2_global, sp = "Mobile MDS2", save = FALSE,
                              y_lab = "mu Mobile PCoA2 (20%)")
plot_count_MDS3_glo <- plot_2(pp = count_MDS3_global, sp = "Mobile MDS3", save = FALSE,
                              y_lab = "mu Mobile PCoA3 (14%)")

plot_cover_MDS1_glo <- plot_2(pp = cover_MDS1_global, sp = "Sessile MDS1", save = FALSE,
                              y_lab = "mu Sessile PCoA1 (19%)")
plot_cover_MDS2_glo <- plot_2(pp = cover_MDS2_global, sp = "Sessile MDS2", save = FALSE,
                              y_lab = "mu Sessile PCoA2 (16%)")
plot_cover_MDS3_glo <- plot_2(pp = cover_MDS3_global, sp = "Sessile MDS3", save = FALSE,
                              y_lab = "mu Sessile PCoA3 (13%)")

plot_count_MDS1_site <- plot_2(pp = site_count_MDS1, sp = "Mobile MDS1", save = TRUE,
                               filename = "supplement/MV_TS_site/count_MDS1_site")
plot_count_MDS2_site <- plot_2(pp = site_count_MDS2, sp = "Mobile MDS2", save = TRUE,
                               filename = "supplement/MV_TS_site/count_MDS2_site")
plot_count_MDS3_site <- plot_2(pp = site_count_MDS3, sp = "Mobile MDS3", save = TRUE,
                               filename = "supplement/MV_TS_site/count_MDS3_site")

plot_cover_MDS1_site <- plot_2(pp = site_cover_MDS1, sp = "Sessile MDS1", save = TRUE,
                               filename = "supplement/MV_TS_site/cover_MDS1_site")
plot_cover_MDS2_site <- plot_2(pp = site_cover_MDS2, sp = "Sessile MDS2", save = TRUE,
                               filename = "supplement/MV_TS_site/cover_MDS2_site")
plot_cover_MDS3_site <- plot_2(pp = site_cover_MDS3, sp = "Sessile MDS3", save = TRUE,
                               filename = "supplement/MV_TS_site/cover_MDS3_site")


#### Add small lines at the bottom everytime deriv2 != 0
#### the times at which deriv2 != 0 are identified from the
#### deriv2 figures in mat sup


plot_count_MDS1_glo <- plot_count_MDS1_glo +
  geom_segment(
    data = data.frame(x = c(1983,1989,1992,1997,2001,2005,2010,2014,2019)),
    aes(x = x, xend = x, y = -Inf, yend = -0.3),
    color = "black",
    inherit.aes = FALSE,
    linewidth = 0.75,
    alpha = 0.5
  )


plot_count_MDS3_glo <- plot_count_MDS3_glo +
  geom_segment(
    data = data.frame(x = 1988),
    aes(x = x, xend = x, y = -Inf, yend = -0.3),
    color = "black",
    inherit.aes = FALSE,
    linewidth = 0.75,
    alpha = 0.5
  )



plot_cover_MDS1_glo <- plot_cover_MDS1_glo +
  geom_segment(
    data = data.frame(x = c(1988,2000)),
    aes(x = x, xend = x, y = -Inf, yend = -0.295),
    color = "black",
    inherit.aes = FALSE,
    linewidth = 0.75,
    alpha = 0.5
  )



#### combine some of the figures for the manuscript

fig3 <- ggarrange(plot_count_MDS1_glo,plot_count_MDS2_glo,plot_count_MDS3_glo,
                  labels = c("a)", "b)", "c)"),
                  #hjust = -0.2, vjust = -0.075,
                  font.label = list(size = 15, face = "bold.italic"),
                  nrow = 1, ncol = 3) %>% 
  annotate_figure(top = text_grob("", color = "white", size = 20), 
  )
fig3

fig5 <- ggarrange(plot_cover_MDS1_glo,plot_cover_MDS2_glo,plot_cover_MDS3_glo,
                  labels = c("a)", "b)", "c)"),
                  #hjust = -0.2, vjust = -0.075,
                  font.label = list(size = 15, face = "bold.italic"),
                  nrow = 1, ncol = 3) %>% 
  annotate_figure(top = text_grob("", color = "white", size = 20), 
  )
fig5

ggsave("time_series_outputs/figures/count_global.tiff",
       plot = fig3, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")

#### combine some of the figures for the manuscript ####################################

### make arrows count ###
dum_data2 <- as.data.frame(cbind(x=c(1,1), y=c(1,8.5)))

arrow_count_MDS1 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "#00A08A")+
  geom_textpath(aes(label = "Mobile PCoA1 (52% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "Empty"), 
            fontface = "italic", size = 4, colour = "#00A08A")+
  geom_text(aes(x = 1, y = 9.25, label = "Littorina \n Nucella"),
            fontface = "italic", size = 4, colour = "#00A08A")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.4))+
  theme_void()+theme(legend.position="none")

arrow_count_MDS2 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "#00A08A")+
  geom_textpath(aes(label = "Mobile PCoA2 (20% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "Amphipoda \n Nucella"),
            fontface = "italic", size = 4, colour = "#00A08A")+
  geom_text(aes(x = 1, y = 9.25, label = "Littorina"),
            fontface = "italic", size = 4, colour = "#00A08A")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.4))+
  theme_void()+theme(legend.position="none")

arrow_count_MDS3 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "#00A08A")+
  geom_textpath(aes(label = "Mobile PCoA3 (14% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "Amphipoda"), 
            fontface = "italic", size = 4, colour = "#00A08A")+
  geom_text(aes(x = 1, y = 9.25, label = "Nucella \n C. maenas"),
            fontface = "italic", size =4, colour = "#00A08A")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.4))+
  theme_void()+theme(legend.position="none")


fig3 <- ggarrange(
  arrow_count_MDS1,
  plot_count_MDS1_glo,
  arrow_count_MDS2,
  plot_count_MDS2_glo,
  arrow_count_MDS3,
  plot_count_MDS3_glo,
  labels = c("", "a)", "", "b)", "", "c)"),
  label.x = -0.05,
  #hjust = -0.2, vjust = -0.075,
  font.label = list(size = 15, face = "bold.italic"),
  nrow = 1,
  ncol = 6,
  widths = c(1,4,1,4,1,4)
) %>%
  annotate_figure(top = text_grob("", color = "white", size = 20),)
fig3


ggsave("time_series_outputs/figures/count_global_arrow.tiff",
       plot = fig3, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")





### Cover ###


arrow_cover_MDS1 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "#046C9A")+
  geom_textpath(aes(label = "Mobile PCoA1 (19% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "Arthopyrenia \n Cyanobacteria"), 
            fontface = "italic", size = 4, colour = "#046C9A")+
  geom_text(aes(x = 1, y = 9.25, label = "S. balanoides \n M. edulis \n Ascophyllum \n fucus"),
            fontface = "italic", size = 4, colour = "#046C9A")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.5))+
  theme_void()+theme(legend.position="none")

arrow_cover_MDS2 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "#046C9A")+
  geom_textpath(aes(label = "Mobile PCoA2 (16% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "C. crispus \n Mastocarpus"), 
            fontface = "italic", size = 4, colour = "#046C9A")+
  geom_text(aes(x = 1, y = 9.25, label = "S. balanoides \n Arthopyrenia"), 
            fontface = "italic", size = 4, colour = "#046C9A")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.5))+
  theme_void()+theme(legend.position="none")

arrow_cover_MDS3 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "#046C9A")+
  geom_textpath(aes(label = "Mobile PCoA3 (13% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "Ascophyllum \n Mastocarpus"), 
            fontface = "italic", size = 4, colour = "#046C9A")+
  geom_text(aes(x = 1, y = 9.25, label = "M. edulis"), 
            fontface = "italic", size =4, colour = "#046C9A")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.5))+
  theme_void()+theme(legend.position="none")




fig5 <- ggarrange(
  arrow_cover_MDS1,
  plot_cover_MDS1_glo,
  arrow_cover_MDS2,
  plot_cover_MDS2_glo,
  arrow_cover_MDS3,
  plot_cover_MDS3_glo,
  labels = c("", "a)", "", "b)", "", "c)"),
  #hjust = -0.2, vjust = -0.075,
  label.x = -0.05, 
  font.label = list(size = 15, face = "bold.italic"),
  nrow = 1,
  ncol = 6,
  widths = c(1.2,3.9,1.2,3.9,1.1,3.9)
) %>%
  annotate_figure(top = text_grob("", color = "white", size = 20),)
fig5


ggsave("time_series_outputs/figures/cover_global_arrow.tiff",
       plot = fig5, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")


### NOAA PCoA  ##########################################################
#____________________________________________________________________________

# Plot arrow


lab_arrows_noaa <- c("Subtidal PCoA1 (23% var)","Subtidal PCoA1 (23% var)",
                     "Subtidal PCoA2 (11% var)", "Subtidal PCoA2 (11% var)") 

fleche_noaa_PCoA <- ggplot(dum_data, aes(x = x, y = y, color = axis)) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 3)+
  geom_textpath(aes(label = lab_arrows_noaa), vjust = -0.5, size = 7)+
  scale_color_manual(values = c("darkgoldenrod","darkred"))+
  geom_text(aes(x = 1, y = 0.5, label = neg_cor_spe_subtidal_PCoA1), size = 4, colour = "darkgoldenrod")+
  geom_text(aes(x = 1, y = 9.25, label = pos_cor_spe_subtidal_PCoA1), size = 4, colour = "darkgoldenrod")+
  geom_text(aes(x = 2, y = 9.2, label = pos_cor_spe_subtidal_PCoA2), size = 4, colour = "darkred")+
  geom_text(aes(x = 2, y = 0.5, label = neg_cor_spe_subtidal_PCoA2), size = 4, colour = "darkred")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.5))+
  theme_void()+theme(legend.position="none")
  
arrow_noaa_MDS1 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "darkgoldenrod")+
  geom_textpath(aes(label = "Subtidal PCoA1 (23% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 0.5, label = "S.scombrus \n A.lupus \n M.americanus"),
            fontface = "italic", size = 4, colour = "darkgoldenrod")+
  geom_text(aes(x = 1, y = 9.25, label = "M.weitzmani \n Sepiolidae \n Shrimps"),
            fontface = "italic", size =4, colour = "darkgoldenrod")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.5))+
  theme_void()+
  theme(legend.position="none")

arrow_noaa_MDS2 <- ggplot(dum_data2, aes(x = c(1,1), y = c(1,8.5))) + 
  geom_path(arrow = arrow(ends = "both", type = "open", angle = 25, length = unit(0.5, "cm")),
            size = 2, color = "darkred")+
  geom_textpath(aes(label = "Subtidal PCoA2 (11% var)"), vjust = -0.5, size = 5)+
  geom_text(aes(x = 1, y = 9.25, label = "Shrimp \n Octopoda"), 
            fontface = "italic", size = 4, colour = "darkred")+
  geom_text(aes(x = 1, y = 0.5, label = "Large fish"),
            fontface = "italic", size =4, colour = "darkred")+
  scale_y_continuous(limits = c(-1, 10))+
  scale_x_continuous(limits = c(0.5, 1.5))+
  theme_void()+theme(legend.position="none")



# plot TS with noaa PCoA 1 et 2

p1 <- ggplot() +
  #geom_line(data = noaa_MDS2, aes(y = estimate, x = year), color = "darkred", size = 1)+
  #geom_ribbon(data = noaa_MDS2, aes(y = estimate, x = year, ymin = lower_CI, ymax = upper_CI), fill = "darkred", color = "darkred", alpha = 0.15)+
  geom_line(data = noaa_MDS1, aes(y =estimate__, x = year), color = "goldenrod", size = 1)+
  geom_ribbon(data = noaa_MDS1, aes(y = estimate__, x = year, ymin = lower__, ymax = upper__),
              fill = "goldenrod", color = "goldenrod", alpha = 0.25)+
  # ylab("")+
  ylab("Subtidal PCoA1 (37% var)")+
  xlab("Time")+
  # geom_vline(xintercept = 2000.7, color ="red", linetype = "longdash", size = 1.5) +
  # geom_rect(aes(xmin = 1999.4, xmax = 2002, ymin=-Inf, ymax=Inf),
  #           fill = "red", alpha = 0.1)+
  ### figure for noaa-keen validation
  geom_vline(xintercept = 2010, color ="red", linetype = "longdash", size = 1.5) +
  geom_rect(aes(xmin = 2009, xmax = 2011, ymin=-Inf, ymax=Inf),
            fill = "red", alpha = 0.1)+
  theme_bw()+
  theme(panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=11),
        axis.title=element_text(size=13),
        #axis.title.y = element_blank()
  )
p1


p2 <- ggplot() +
  geom_line(data = noaa_MDS2, aes(y = estimate__, x = year), color = "darkred", size = 1)+
  geom_ribbon(data = noaa_MDS2, aes(y = estimate__, x = year, ymin = lower__, ymax = upper__),
              fill = "darkred", color = "darkred", alpha = 0.25)+
  #geom_line(data = noaa_MDS1, aes(y =estimate, x = year), color = "goldenrod", size = 1)+
  #geom_ribbon(data = noaa_MDS1, aes(y = estimate, x = year, ymin = lower_CI, ymax = upper_CI), fill = "goldenrod", color = "goldenrod", alpha = 0.15)+
  # ylab("")+
  ylab("Subtidal PCoA3 (14% var)")+
  xlab("Time")+
  # geom_vline(xintercept = c(1990.2,2010.2), color ="red", linetype = "longdash", size = 1.5) +
  # geom_rect(aes(xmin = c(1987,2009.5), xmax = c(1992.7,2011.7), ymin=-Inf, ymax=Inf), 
  #           fill = "red", alpha = 0.1)+
  geom_vline(xintercept = c(1990.2,2000), color ="red", linetype = "longdash", size = 1.5) +
  geom_rect(aes(xmin = c(1987,1999), xmax = c(1992.7,2002), ymin=-Inf, ymax=Inf), 
            fill = "red", alpha = 0.1)+
  theme_bw()+
  theme(panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=11),
        axis.title=element_text(size=13),
        #axis.title.y = element_blank()
  )
p2



# put them together


noaa_PCoA_1_2_fl <- ggarrange(arrow_noaa_MDS1, p1, arrow_noaa_MDS2, p2, 
                              labels = c("", "a)", "", "b)"),
                              #hjust = -0.2, vjust = -0.075,
                              font.label = list(size = 15, face = "bold.italic"),
                              nrow = 1,
                              ncol = 4,
                              widths = c(1,5,1,5)
)

# noaa_PCoA_1_2_fl <- ggarrange(p1, 
#                               p2, 
#                               ncol = 2, widths = c(5,5))
# 
#  ggsave("time_series_outputs/figures/NOAA_TS_PCoA-KEEN.tiff",
#         plot = noaa_PCoA_1_2_fl, device = "tiff", scale = 2,
#         width = 1920, height = 991, bg = "white", units = "px", dpi = 300,
#         compression = "lzw")


ggsave("time_series_outputs/figures/NOAA_TS_PCoA1-2_arrow.tiff",
       plot = noaa_PCoA_1_2_fl, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")



