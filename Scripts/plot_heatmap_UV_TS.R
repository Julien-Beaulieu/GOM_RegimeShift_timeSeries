# This script is to plot a heat map with the synchrony between changes in 
# intertidal species abundances and subtidal regime shifts. The synchrony was
# determined based on the codels and graphs in the script plot_UV_TS.R
# Contributor: Julien Beaulieu

# Load packages
library(ggplot2)
library(tidyverse)
library(data.table)
library(reshape2)
library(geomtextpath)
library(ggpubr)

# 1. Load and format the 2nd derivative values
deriv_values <- read.csv("time_series_outputs/values_deriv2_2.csv")
deriv <- deriv_values %>% select(-X) %>% 
  mutate(specie = 
           ifelse(str_detect(var,"Ascophyllum"),"A. nodusum",
                  ifelse(str_detect(var,"maenas"),"C. maenas",
                         ifelse(str_detect(var,"Chondrus"),"C. crispus",
                                ifelse(str_detect(var,"Fucus"),"Fucus sp.",
                                       ifelse(str_detect(var,"Littorina"),"Littorina sp.",
                                              ifelse(str_detect(var,"Mastocarpus"),"M. stellatus",
                                                     ifelse(str_detect(var,"ytilus"),"M. edulis",
                                                            ifelse(str_detect(var,"Nucella"),"Nucella sp.",
                                                                   ifelse(str_detect(var,"Semibalanus"),"S. balanoides",
                                                                          ifelse(str_detect(var,"Amphipoda"),"Amphipoda",
                                                                                 ifelse(str_detect(var,"Isopoda"),"Isopoda","-")))))))))))) %>% 
  mutate(site = ifelse(str_detect(var,"NW"),"NW",
                       ifelse(str_detect(var,"SW"),"SW",
                              ifelse(str_detect(var,"NE"),"NE",
                                     ifelse(str_detect(var,"SE"),"SE","Babb's Cove")))))

# 2. Create data frame automatically (Matches names in 'deriv')
species_list <- c("Amphipoda","Isopoda","C. maenas","Littorina sp.",
                  "Nucella sp.", "M. edulis", "S. balanoides", "C. crispus",
                  "M. stellatus", "A. nodusum","Fucus sp.")
sites_list <- c("Babb's Cove","SE","NE","SW","NW")

data_heatmap <- expand.grid(specie = species_list, site = sites_list) %>%
  mutate(nb = 0)

# 3. Loop to count regime shifts
for (i in 1:nrow(data_heatmap)) {
  temp <- deriv %>% filter(specie == data_heatmap$specie[i] &
                             site == data_heatmap$site[i])
  if(nrow(temp) > 0){
    for (j in 1:min(3, nrow(temp))) {
      if((temp$upper[j] * temp$lower[j] > 0) |
         (temp$upper_zi[j] * temp$lower_zi[j] > 0)){
        data_heatmap$nb[i] = data_heatmap$nb[i] + 1
      }
    }
  }
}

# 4. Set Factors for Plotting Order
data_heatmap$specie <- factor(data_heatmap$specie, levels = species_list)
data_heatmap$site <- factor(data_heatmap$site, levels = sites_list)

# 5. Plot Main Heatmap (Synchrony Count)
heatmap <- ggplot(data = data_heatmap, aes(x = specie, y = site, fill = as.factor(nb)))+
  geom_tile(color = "black", linewidth = 0.2)+
  scale_fill_brewer(palette = "Greys", name = "Synchrony") +
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14))

# 6. Function for Categorical Directional Heatmaps
heatmap_plot <- function(der = deriv_RS1, dat = data_heatmap){
  dat$nb = 0
  for (i in 1:nrow(dat)) {
    temp <- der %>% filter(specie == dat$specie[i] & site == dat$site[i])
    if(nrow(temp) > 0){
      if(temp$upper * temp$lower > 0){
        if(temp$upper < 0){dat$nb[i] = -1}
        if(temp$lower > 0){dat$nb[i] = 1}
      }
      # Logic for Zero-Inflation discrepancies
      if(temp$upper_zi * temp$lower_zi > 0){
        if((temp$lower_zi > 0 & dat$nb[i] == 1) | (temp$upper_zi < 0 & dat$nb[i] == -1)){
          dat$nb[i] = 0
          print(paste("WARNING discrepancy at", dat$specie[i], dat$site[i]))
        } else if(temp$upper_zi < 0){
          dat$nb[i] = 1
        } else if(temp$lower_zi > 0){
          dat$nb[i] = -1
        }
      }
    }
  }
  
  dat$specie <- factor(dat$specie, levels = species_list)
  dat$site <- factor(dat$site, levels = sites_list)
  
  ggplot(data = dat, aes(x = specie, y = site, fill = as.factor(nb)))+
    geom_tile(color = "gray80")+
    scale_fill_manual(
      name = "Offset synchrony",
      values = c("-1" = "darkblue", "0" = "white", "1" = "red"),
      labels = c("-1" = "Lower", "0" = "Stable", "1" = "Higher"),
      drop = FALSE
    ) +
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

### Make 3 heatmaps (1/RS) with color code for + and - effects #################

deriv_RS1 <- deriv %>% filter(regime_shift == 1)
deriv_RS2 <- deriv %>% filter(regime_shift == 2)
deriv_RS3 <- deriv %>% filter(regime_shift == 3)

#### run the plots #######################

plot_RS1 <- heatmap_plot()
plot_RS2 <- heatmap_plot(der = deriv_RS2)
plot_RS3 <- heatmap_plot(der = deriv_RS3)

heatmap_plot2 <- ggarrange(plot_RS1,plot_RS2,plot_RS3,
                           labels = c("a) 1990 Regime Shift",
                                      "b) 2000 Regime Shift",
                                      "c) 2010 Regime Shift"),
                           ncol = 3, nrow = 1,
                           font.label = list(size = 13, face = "bold.italic"),
                           common.legend = TRUE, legend = "bottom",
                           vjust = -0.05)
heatmap_plot2_2 <- annotate_figure(heatmap_plot2, top = text_grob("", color = "white", size = 20), 
)

heatmap_plot2_2

ggsave("time_series_outputs/figures/heatmap_V2.tiff",
       plot = heatmap_plot2_2, device = "tiff", scale = 2, 
       width = 1920, height = 991, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")
