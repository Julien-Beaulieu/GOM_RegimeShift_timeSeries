# This script is to produce a time line with the principal perturbation events
# in the GOM from 1980 to 2020. This script also makes conceptual figures of the 
# derivative.
# Contributor: Julien Beaulieu

#load packages

library(brms)
library(ggplot2)
library(grid)
library(ggpubr)
library(janitor)
library(data.table)
library(boot)
library(tidyr)
library(tidyverse)
library(wesanderson)

## Plot time line of Natural history of the GOM ######################################
#___________________________________________________________________________________

noaa <- read.csv("data/NOAA_trawl_data.csv") %>%
  filter(scaled_med_individ_PUE != "NaN") %>% #there are some NaNs in the data but don't include any of our species of interest
  filter(scaled_med_tot_mass_PUE != "NaN") %>%
  mutate(scaled_med_individ_PUE = as.numeric(scaled_med_individ_PUE)) %>%
  mutate(scaled_med_tot_mass_PUE = as.numeric(scaled_med_tot_mass_PUE))

#-Convert to wide dataframe-#
noaa_wide <- dcast(noaa, YEAR ~ SCINAME, value.var = "scaled_med_tot_mass_PUE", mean) %>%
  as_tibble() %>%
  clean_names() %>%
  subset(year > 1981)

#-NaNs are years where no individuals of that species was detected, will replace with 0-#

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
noaa_wide[is.nan(noaa_wide)] <- 0
noaa_wide$year <- as.numeric(noaa_wide$year)

ggplot(data = noaa_wide, aes(x = year, y = gadus_morhua))+
  geom_point() + geom_smooth()
ggplot(data = noaa_wide, aes(x = year, y = homarus_americanus))+
  geom_point() + geom_smooth()
ggplot(data = noaa_wide, aes(x = year, y = brosme_brosme))+
  geom_point() + geom_smooth()
ggplot(data = noaa_wide, aes(x = year, y = lumpenus_lumpretaeformis))+
  geom_point() + geom_smooth()
ggplot(data = noaa_wide, aes(x = year, y = lithodes_maja))+
  geom_point() + geom_smooth()

#Create data with time interval for major events

col_name <- c("Event","Start","End","y","type")
cod <- c("Fisheries collapse","1980","1990",1,"collapse") # ref = noaa data
urchins1 <- c("Rise of urchins","1980","1995",2,"Rise") #steneck
urchins2 <- c("Collapse of urchins", "1995","2000",2,"collapse") #Steneck
reds <- c("Rise of filamentous red algae","2010","2015",4,"Rise") #Jen's paper
lobster <- c("Rise of mesopredators","1997","2020",3,"Rise") #noaa data; start = 2005;; Steneck
SST <- c("Climate change impacts","2010","2020",1,"Climate") # includes roughly temperature and anomalies; from Pershing et al. 2021, Bricknell et al. 2021
history_data <- as.data.frame(rbind(cod,urchins1,urchins2,reds,lobster,SST))
colnames(history_data) <- col_name
history_data$Start <- as.numeric(history_data$Start)
history_data$End <- as.numeric(history_data$End)

# Plot

timeline <- ggplot(data = history_data, aes(x = Start, y = y))+
  scale_color_manual(values = c("darkgreen","darkred","darkblue"))+
  geom_segment(aes(xend = End, yend = y, color = type), size = 5, alpha = 0.6)+
  geom_text(aes(label = Event, vjust = -1, hjust = 0), size = 6)+
  xlab("")+
  ylab("  ")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        text = element_text(size = 15),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.text=element_text(size=22),
  )
timeline

ggsave("time_series_outputs/figures/timeline_GOM.tiff",
       plot = timeline, device = "tiff", scale = 2, 
       width = 1920, height = 500, bg = "white", units = "px", dpi = 300, 
       compression = "lzw")


#### Make figure for hypothesis ###################################

inv_logit <- function(x) 1 / (1 + exp(-x))

sub <- as.data.frame(x = c(1:3000)) %>% 
  rename("x" = `c(1:3000)`) %>% 
  mutate(y = inv_logit(-exp(1.8)*(log(x)-6.75))) %>% 
  ggplot(aes(x=x,y=y))+
  geom_rect(aes(xmin = exp(6.75)-330, xmax = exp(6.75)+430, ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.03)+
  geom_line(linewidth = 2, color = "#972D15")+
  xlab("Time")+
  ylab("Subtidal community")+
  geom_vline(xintercept = exp(6.75), linetype = "dashed", color = "red",
             linewidth = 1.5)+
  geom_vline(xintercept = exp(6.75)-350, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_vline(xintercept = exp(6.75)+450, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_text(aes(x = exp(6.75)-425, y = 0.5, label = "RS start"), angle = 90,
            size = 15)+
  geom_text(aes(x = exp(6.75)+375, y = 0.5, label = "RS end"), angle = 90,
            size = 15)+
  theme_bw()+
  theme(axis.text = element_blank(),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title=element_text(size=40))

sub



inter <- as.data.frame(x = c(1:3000)) %>% 
  rename("x" = `c(1:3000)`) %>% 
  mutate(y = inv_logit(-exp(1.8)*(log(x)-6.75))) %>%
  mutate(y2 = inv_logit(exp(1.8)*(log(x)-7.3))) %>% 
  ggplot(aes(x=x,y=y))+
  geom_rect(aes(xmin = exp(6.75)-330, xmax = exp(6.75)+430, ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.03)+
  geom_line(linewidth = 2, color = "darkgrey", alpha = 0.5)+
  xlab("Time")+
  ylab("Intertidal community")+
  geom_vline(xintercept = exp(6.75), linetype = "dashed", color = "black",
             linewidth = 1.5, alpha = 0.3)+
  geom_vline(xintercept = exp(6.75)-350, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_vline(xintercept = exp(6.75)+450, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  # geom_text(aes(x = exp(6.75)-425, y = 0.5, label = "RS start"), 
  #           size = 15, angle = 90, color = "grey")+
  # geom_text(aes(x = exp(6.75)+375, y = 0.5, label = "RS end"),
  #           size = 15, angle = 90, color = "grey")+
  geom_line(aes(y = y2), linewidth = 2.5, color = "#046C9A")+
  theme_bw()+
  theme(axis.text = element_blank(),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title=element_text(size=40)
        )

inter


inter2 <- as.data.frame(x = c(1:3000)) %>% 
  rename("x" = `c(1:3000)`) %>% 
  mutate(y = inv_logit(-exp(1.8)*(log(x)-6.75))) %>%
  mutate(y2 = inv_logit(exp(1.8)*(log(x)-7.8))) %>% 
  ggplot(aes(x=x,y=y))+
  geom_rect(aes(xmin = exp(6.75)-330, xmax = exp(6.75)+430, ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.03)+
  geom_line(linewidth = 2, color = "darkgrey", alpha = 0.5)+
  xlab("Time")+
  ylab("Intertidal community")+
  geom_vline(xintercept = exp(6.75), linetype = "dashed", color = "black",
             linewidth = 1.5, alpha = 0.3)+
  geom_vline(xintercept = exp(6.75)-350, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_vline(xintercept = exp(6.75)+450, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  # geom_text(aes(x = exp(6.75)-425, y = 0.5, label = "RS start"), 
  #           size = 15, angle = 90, color = "grey")+
  # geom_text(aes(x = exp(6.75)+375, y = 0.5, label = "RS end"),
  #           size = 15, angle = 90, color = "grey")+
  geom_line(aes(y = y2), linewidth = 2.5, color = "#046C9A")+
  theme_bw()+
  theme(axis.text = element_blank(),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title=element_text(size=40)
  )

inter2


ggsave("time_series_outputs/figures/hyp_subtidal.tiff",
       plot = sub, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
       )

ggsave("time_series_outputs/figures/hyp_intertidal.tiff",
       plot = inter, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
       )

ggsave("time_series_outputs/figures/hyp2_intertidal.tiff",
       plot = inter2, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)

ggsave("time_series_outputs/figures/timeline_subtidal.tiff",
       plot = timeline, device = "tiff", scale = 2, 
       width = 1920, height = 591, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)


### plot deriv for second panel ##########################

### calculer driv1 et deriv2

dat <- as.data.frame(x = c(1:3000)) %>%
  rename("x" = `c(1:3000)`) %>%
  mutate(y = inv_logit(-exp(1.8) * (log(x) - 6.75)))

deriv1 <- function(data = dat){
  y <- data %>% select(y)
  y1 <- y[-1,]
  y0 <- y[-nrow(y),]
  
  num = y1-y0
  
  x <- data %>% select(x)
  x1 <- x[-1,]
  x0 <- x[-nrow(x),]
  
  deno = x1-x0
  
  deriv <- as.data.frame(num/deno) %>% mutate(x = x1) %>% 
    rename(y_der = `num/deno`)
  
  return(deriv)

}
  
der1_dat <- deriv1(data=dat)  
der2_dat <- der1_dat %>% rename("y" = y_der) %>% deriv1()  

### plot

der0_plot <- dat %>% filter(x<2010) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line(linewidth = 2)+
  geom_vline(xintercept = exp(6.75)-40, linetype = "dashed", color = "black",
             linewidth = 1.5, alpha = 0.3)+
  geom_hline(yintercept = 0, color = "red")+
  ylab("Community")+
  xlab("Time")+
  geom_vline(xintercept = exp(6.75)-225, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_vline(xintercept = exp(6.75)+125, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  theme_bw()+
  theme(axis.text = element_blank(),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title=element_text(size=50)
  )

der0_plot


der1_plot <- der1_dat %>% filter(x<2010) %>% 
  ggplot(aes(x = x, y = y_der)) +
  geom_line(linewidth = 2)+
  geom_vline(xintercept = exp(6.75)-40, linetype = "dashed", color = "black",
             linewidth = 1.5, alpha = 0.3)+
  geom_hline(yintercept = 0, color = "red", alpha = 1)+
  ylab("First derivative")+
  xlab("Time")+
  geom_vline(xintercept = exp(6.75)-225, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_vline(xintercept = exp(6.75)+125, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_text(aes(x = exp(6.75)-100, y = -0.0008, label = "Fastest change"),
            size = 12, angle = 90, color = "black", alpha = 0.8)+
  theme_bw()+
  theme(axis.text = element_blank(),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title=element_text(size=50)
  )
  
der1_plot


der2_plot <- der2_dat %>% filter(x<2010) %>% 
  ggplot(aes(x = x, y = y_der)) +
  geom_line(linewidth = 2)+
  geom_vline(xintercept = exp(6.75)-40, linetype = "dashed", color = "black",
             linewidth = 1.5, alpha = 0.3)+
  geom_hline(yintercept = 0, color = "red", alpha = 1)+
  ylab("Second derivative")+
  xlab("Time")+
  geom_vline(xintercept = exp(6.75)-225, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_vline(xintercept = exp(6.75)+125, linetype = "dashed", color = "darkgrey",
             linewidth = 1)+
  geom_text(aes(x = exp(6.75)-300, y = 0, label = "Shift in trend (-)"), angle = 90,
            size = 12)+
  geom_text(aes(x = exp(6.75)+50, y = 0, label = "Shift in trend (+)"), angle = 90,
            size = 12)+
  theme_bw()+
  theme(axis.text = element_blank(),
        panel.background = element_blank(), 
        #panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title=element_text(size=50)
  )

der2_plot


## save
ggsave("time_series_outputs/figures/der0_conceptual_plot.tiff",
       plot = der0_plot, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)

ggsave("time_series_outputs/figures/der1_conceptual_plot.tiff",
       plot = der1_plot, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)

ggsave("time_series_outputs/figures/der2_conceptual_plo1.tiff",
       plot = der2_plot, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)


  