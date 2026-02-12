# this script is to plot morinations. The score or the PCoA are from the script 
# multivariate_time_series. The PCoA scores (input) are produced by the script 
# multivariate_time_series
# Contributor: Julien Beaulieu

# Load packages
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)

#load data 
score.noaa <- read.csv("time_series_outputs/ordination_output/score_noaa_pcoa_site.csv")
spp.score <- read.csv( "time_series_outputs/ordination_output/score_noaa_specie.csv")
spp.score <- spp.score[,-1]

cover_merge <- read.csv("time_series_outputs/ordination_output/score_PCoAsite_cover.csv")
spp.score.cover <- read.csv("time_series_outputs/ordination_output/spp_score_cover.csv")
spp.score.cover <- spp.score.cover[,-1]

count_merge <- read.csv("time_series_outputs/ordination_output/score_PCoAsite_count.csv")
spp.score.count <- read.csv("time_series_outputs/ordination_output/spp_score_count.csv")
spp.score.count <- spp.score.count[,-1]


# PLOT PCoA #####################################################################
###_______________________________________________________________________________
score.noaa$year <- as.numeric(score.noaa$year)

score.noaa$period <- NA

for (i in c(1:nrow(score.noaa))) {
  if(score.noaa$year[i] <= 1990){score.noaa$period[i] = "1"}
  if(score.noaa$year[i] > 1990){score.noaa$period[i] = "2"}
  if(score.noaa$year[i] > 2000){score.noaa$period[i] = "3"}
  if(score.noaa$year[i] > 2010){score.noaa$period[i] = "4"}
}

score.noaa$period <- as.factor(score.noaa$period)
col <- c("blue","purple","pink","red")

spp.score$no <- as.character(c(1:nrow(spp.score)))

plot_pcoa_noaa <- ggplot()+ 
  geom_point(data = score.noaa, aes(y = MDS2, x = MDS1, fill = year), shape = 21, alpha = 0.7, size = 5)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  #stat_ellipse(data = score.noaa, geom="path", aes(y = MDS2, x = MDS1, group = period), #?lipses pointill
  #             linetype = "dotted",
  #             alpha = 1, 
  #             show.legend = FALSE, 
  #             level = 0.95) +
  #stat_ellipse(data = score.noaa, geom="polygon", aes(y = MDS2, x = MDS1, fill = period), #?lipses
  #             alpha = 0.12, 
  #             show.legend = FALSE, 
  #             level = 0.95) +
  #scale_fill_manual(values = col)+
  geom_text_repel(data = spp.score, aes(x = MDS1, y = MDS2, label = species), max.overlaps = 25, point.padding = 0.1, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (23%)") + 
  ylab("PCoA 2 (11%)") +
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text = element_text(size=16),
        legend.title = element_text(size=17),
        #legend.position = "top",
        #legend.direction =  "vertical"
  )
plot_pcoa_noaa


#save PCoA plot
#ggsave("time_series_outputs/figures/supplement/PCoA_NOAA_mat_sup2.tiff", plot = plot_pcoa_noaa, device = "tiff", scale = 2, width = 1920, height = 991, bg = "white",
#       units = "px", dpi = 300)

#save species list and numbers as a table

spp_names <- as.data.frame(spp.score[,c(1,5)])
colnames(spp_names) <- c("Taxa","Numerical ID") 

#save species list to put in supplement material -- Will have to copy paste because R suck at saving tables as images
#write.csv(spp_names, file = "time_series_outputs/figures/supplement/list_spp_noaa_pcoa.csv",
#          row.names = F)


### PCOA COVER #################################################################

plot_pcoa_cover <- ggplot()+ 
  geom_point(data = cover_merge, aes(y = MDS2, x = MDS1, fill = YEAR, shape = SITE), shape = 21, alpha = 0.7, size = 3)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score.cover, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score.cover, aes(x = MDS1, y = MDS2, label = species), max.overlaps = 25, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (19%)") + 
  ylab("PCoA 2 (16%)") +
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text = element_text(size=16),
        legend.title = element_text(size=17),
        #legend.position = "top",
        #legend.direction =  "vertical"
  )
plot_pcoa_cover

plot_pcoa_cover2 <- ggplot()+ 
  geom_point(data = cover_merge, aes(y = MDS3, x = MDS1, fill = YEAR, shape = SITE), shape = 21, alpha = 0.7, size = 3)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score.cover, aes(x = 0, xend = MDS1, y = 0, yend = MDS3), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score.cover, aes(x = MDS1, y = MDS3, label = species), max.overlaps = 25, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (19%)") + 
  ylab("PCoA 3 (13%)") +
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text = element_text(size=16),
        legend.title = element_text(size=17),
        #legend.position = "top",
        #legend.direction =  "vertical"
  )
plot_pcoa_cover2

ggsave("time_series_outputs/figures/supplement/PCoA_cover_mat_sup.tiff", plot = plot_pcoa_cover, device = "tiff", scale = 2, width = 1920, height = 991, bg = "white",
       units = "px", dpi = 300)

### PCOA COUNT ########################################################

plot_pcoa_count <- ggplot()+ 
  geom_point(data = count_merge, aes(y = MDS2, x = MDS1, fill = YEAR, shape = SITE), shape = 21, alpha = 0.7, size = 4)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score.count, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score.count, aes(x = MDS1, y = MDS2, label = species), max.overlaps = 25, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (52%)") + 
  ylab("PCoA 2 (20%)") +
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text = element_text(size=16),
        legend.title = element_text(size=17),
        #legend.position = "top",
        #legend.direction =  "vertical"
  )
plot_pcoa_count


plot_pcoa_count2 <- ggplot()+ 
  geom_point(data = count_merge, aes(y = MDS3, x = MDS1, fill = YEAR, shape = SITE), shape = 21, alpha = 0.7, size = 4)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score.count, aes(x = 0, xend = MDS1, y = 0, yend = MDS3), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score.count, aes(x = MDS1, y = MDS3, label = species), max.overlaps = 25, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (52%)") + 
  ylab("PCoA 3 (14%)") +
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.text = element_text(size=16),
        legend.title = element_text(size=17),
        #legend.position = "top",
        #legend.direction =  "vertical"
  )
plot_pcoa_count2

ggsave("time_series_outputs/figures/supplement/PCoA_count_mat_sup.tiff", plot = plot_pcoa_count, device = "tiff", scale = 2, width = 1920, height = 991, bg = "white",
       units = "px", dpi = 300)
