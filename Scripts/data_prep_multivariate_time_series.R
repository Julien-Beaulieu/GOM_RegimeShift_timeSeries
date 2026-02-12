# This script is to prepare the data and ordinations to fit multivariate time-series 
# (GAM) on intertidal data and also subtidal data from NOAA and compare the rates of changes.
# Contributor: Julien beaulieu

# 0: load packages
# 1: get data ready
# 1.1: NOAA
# 1.2: intertidal
# 2: data visualisation
# 3: NOAA time serie
# 4: intertidal time serie


## 0: load packages and data ################################################
##________________________________________________________________________

## Packages

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
library(ggrepel)

## 1: Data ##############################################################
#______________________________________________________________________


### 1.1: NOAA


noaa <- read_csv("data/NOAA_trawl_data.csv") %>%
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

# sum species with the same genus
noaa_wide$sum_alosa <- rowSums(noaa_wide[,c("alosa_aestivalis", "alosa_mediocris", "alosa_pseudoharengus","alosa_sapidissima")])
noaa_wide$sum_ammodytes <- rowSums(noaa_wide[,c("ammodytes_americanus","ammodytes_dubius")])
noaa_wide$sum_argentina <- rowSums(noaa_wide[,c("argentina_silus","argentina_striata")])
noaa_wide$sum_cancer <- rowSums(noaa_wide[,c("cancer_borealis","cancer_irroratus")])
noaa_wide$sum_urophycis <- rowSums(noaa_wide[c("urophycis_chesteri", "urophycis_chuss","urophycis_regia",
                                               "urophycis_sp","urophycis_tenuis")])

# remove species that were seen in less than 10 different years
noaa_cut <- noaa_wide
to_rm <- as.data.frame(c(1:ncol(noaa_wide)))
colnames(to_rm)[1] <- "col"
to_rm$n_obs <- 0
for (i in c(2:ncol(noaa_wide))) {
  for (j in c(1:nrow(noaa_wide))){
    if (noaa_wide[j,i] != 0){ 
      to_rm$n_obs[i] <- to_rm$n_obs[i] + 1 
    }
  }
}

for (i in c(2:nrow(to_rm))) {
  if (to_rm$n_obs[i] < 10){
    to_rm$col[i] <- 0 
  }
}
to_rm <- filter(to_rm, col != 0)
to_rm <- to_rm$col
noaa_cut <- noaa_cut[,to_rm]

# remoove species with same genus

drop <-  c("alosa_aestivalis", "alosa_mediocris", "alosa_pseudoharengus",
           "alosa_sapidissima", "ammodytes_americanus","ammodytes_dubius", 
           "argentina_silus","argentina_striata", "cancer_borealis","cancer_irroratus",
           "urophycis_chesteri", "urophycis_chuss","urophycis_regia",
           "urophycis_sp","urophycis_tenuis")
noaa_cut <- noaa_cut[,!(names(noaa_cut) %in% drop)]



### 1.2: intertidal data ####################################################
##______________________________________________________________________

#load
intertidal_data <- read_rds("data/combined_intertidal_abundance.RDS")

# make a unique ID column with site-transect-level-year

intertidal_data <- unite(intertidal_data, "uniqueID", c("SITE","INTERTIDAL_TRANSECT","LEVEL","YEAR"), sep = "-", remove = F)

#remove malaga cut

intertidal_data <- filter(intertidal_data, SITE != "Malaga Cut")

# separate cover and count

intertidal_cover <- filter(intertidal_data, MEASURE == "Percent_Cover")
intertidal_count <- filter(intertidal_data, MEASURE == "Count")

# for each separate in a df for community and a df for variables

cover_community <- intertidal_cover[,c("uniqueID", "ORGANISM","VALUE")]
count_community <- intertidal_count[,c("uniqueID", "ORGANISM","VALUE")]


var <- unique(intertidal_data[,c("uniqueID","SITE","INTERTIDAL_TRANSECT","LEVEL","YEAR")])
var$exposed <- NA
for(i in c(1:nrow(var))){
  if (var$SITE[i] == "Babb's Cove"){
    var$exposed[i] = "1"
  }
  if (var$SITE[i] == "NW Appledore"){
    var$exposed[i] = "0"
  }
  if (var$SITE[i] == "NE Appledore"){
    var$exposed[i] = "1"
  }
  if (var$SITE[i] == "SE Appledore"){
    var$exposed[i] = "1"
  }
  if (var$SITE[i] == "Malaga Cut"){
    var$exposed[i] = "0"
  }
  if (var$SITE[i] == "SW Appledore"){
    var$exposed[i] = "0"
  }
}


#-Convert to wide dataframe-#
cover_wide <- dcast(cover_community, uniqueID ~ ORGANISM, value.var = "VALUE", mean) %>%
  as_tibble() %>%
  clean_names()
count_wide <- dcast(count_community, uniqueID ~ ORGANISM, value.var = "VALUE", mean) %>%
  as_tibble() %>%
  clean_names()

#-NaNs are years where no individuals of that species was detected, will replace with 0-#

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
cover_wide[is.nan(cover_wide)] <- 0

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
count_wide[is.nan(count_wide)] <- 0

# remoove species that were seen in less than 100 times (2.5% of the samples) 
cover_cut <- cover_wide
to_rm <- as.data.frame(c(1:ncol(cover_wide)))
colnames(to_rm)[1] <- "col"
to_rm$n_obs <- 0
for (i in c(2:ncol(cover_wide))) {
  for (j in c(1:nrow(cover_wide))){
    if (cover_wide[j,i] != 0){ 
      to_rm$n_obs[i] <- to_rm$n_obs[i] + 1 
    }
  }
}

for (i in c(2:nrow(to_rm))) {
  if (to_rm$n_obs[i] < 100){
    to_rm$col[i] <- 0 
  }
}

to_rm <- filter(to_rm, col != 0)
to_rm <- to_rm$col
cover_cut <- cover_cut[,to_rm]
cover_wide <- cover_cut

count_cut <- count_wide
to_rm <- as.data.frame(c(1:ncol(count_wide)))
colnames(to_rm)[1] <- "col"
to_rm$n_obs <- 0
for (i in c(2:ncol(count_wide))) {
  for (j in c(1:nrow(count_wide))){
    if (count_wide[j,i] != 0){ 
      to_rm$n_obs[i] <- to_rm$n_obs[i] + 1 
    }
  }
}

for (i in c(2:nrow(to_rm))) {
  if (to_rm$n_obs[i] < 100){
    to_rm$col[i] <- 0 
  }
}

to_rm <- filter(to_rm, col != 0)
to_rm <- to_rm$col
count_cut <- count_cut[,to_rm]
count_wide <- count_cut



### 3: NOAA timeseries #######################################################
##____________________________________________________________________________

# make ordinations (PCoA) using bray-curtis distance ###################################
##______________________________________________________________________________

noaa_cut <- noaa_cut %>% remove_rownames %>% column_to_rownames(var="year") 


# Calculate PCoA

community <- decostand(noaa_cut, method = "hellinger") 



#PCoA
pcoa <- capscale(community~1, distance = "euclidean") # 25 et 19 %
summary(pcoa) #23%, 11%, 7%

#extraire scores sites

score.pcoa.sites <- as.data.frame(scores(pcoa, choices = c(1:3), display = "sites"))
setDT(score.pcoa.sites, keep.rownames = "year")[]
score.noaa <- score.pcoa.sites

#extraire score especes

spp.score <- data.frame(scores(pcoa, display = "species", choices = c(1:3)))
setDT(spp.score, keep.rownames = "species")[]

write.csv(score.pcoa.sites, "time_series_outputs/ordination_output/score_noaa_pcoa_site.csv")
write.csv(spp.score, "time_series_outputs/ordination_output/score_noaa_specie.csv")

### Run model on coordinates 

score.pcoa.sites$year <- as.numeric(score.pcoa.sites$year)
noaa_pcoa_mod <- bf( mvbind(MDS1,MDS2,MDS3) ~ s(year), family = gaussian())
get_prior(noaa_pcoa_mod, data = score.pcoa.sites)

prior_noaa_pcoa_mod <- c(prior(normal(0,1), class = "b", resp = ""),
                         prior(normal(0,1), class = "b", resp ="MDS1"),
                         prior(normal(0,1), class = "b", resp = "MDS2"),
                         prior(normal(0,1), class = "b", resp = "MDS3"),
                         prior(normal(0,1), class = "Intercept", resp = ""),
                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS1"),
                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS2"),
                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS3"),
                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS1"),
                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS2"),
                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS3"),
                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS1"),
                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS2"),
                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS3"))

## fit !!!!

#fit <- brm(noaa_pcoa_mod,
#           prior = prior_noaa_pcoa_mod,
#           data = score.pcoa.sites,
#           warmup = 2000, iter = 7000,
#           cores = 3, chains = 3,
#           control = list(adapt_delta = 0.99))

#saveRDS(fit, file = "time_series_outputs/noaa_pcoa_GS_3.rds")



### Intertidal ordinations timeseries #######################################
###______________________________________________________________________

# sum littorina species because we don't trust ID

count_wide$sum_littorina <- count_wide$littorina_littorea + count_wide$littorina_obtusata + count_wide$littorina_saxatilis
count_wide[c("littorina_littorea","littorina_obtusata","littorina_saxatilis")] <- list(NULL)

# sum fucus because we don't trust ID

cover_wide[c("fucus_distichus","fucus_spiralis", "fucus_vesiculosus", "shell_hash")] <- list(NULL)

## Ordinations and score extraction for count and cover data

count_wide <- count_wide %>% remove_rownames %>% column_to_rownames(var="unique_id") 
count_wide <- decostand(count_wide, method = "hellinger") 

cover_wide <- cover_wide %>% remove_rownames %>% column_to_rownames(var="unique_id") 
cover_wide <- decostand(cover_wide, method = "hellinger") 

# remove anurida_maritima because low confidence in ID

count_wide["anurida_maritima"] <- list(NULL)

#Remove sessile organism from count data

count_wide[c("mytilus_edulis","semibalanus_balanoides","modiolus_modiolus","testudinalia_testudinalis")] <- list(NULL)

# remove bare rocks and black zone from % cover

cover_wide[c("bare_rock","black_zone","little_round_green_things")] <- list(NULL)

#PCoA
pcoa_count <- capscale(count_wide~1, distance = "euclidean") 
summary(pcoa_count) #52, 20, 14


pcoa_cover <- capscale(cover_wide~1, distance = "euclidean") 
summary(pcoa_cover) # 19% ,16% , 13%

#extraire scores sites

score.count <- as.data.frame(scores(pcoa_count, choices = c(1:3), display = "sites"))
setDT(score.count, keep.rownames = "unique_id")[]

score.cover <- as.data.frame(scores(pcoa_cover, choices = c(1:3), display = "sites"))
setDT(score.cover, keep.rownames = "unique_id")[]


#extraire score especes

spp.score.count <- data.frame(scores(pcoa_count, display = "species", choices = c(1:3)))
setDT(spp.score.count, keep.rownames = "species")[]

spp.score.cover <- data.frame(scores(pcoa_cover, display = "species", choices = c(1:3)))
setDT(spp.score.cover, keep.rownames = "species")[]


# merge scores with other variables

count_merge <- unique(merge(score.count, var, by.x = "unique_id", by.y = "uniqueID"))
cover_merge <- unique(merge(score.cover, var, by.x = "unique_id", by.y = "uniqueID"))  

write.csv(cover_merge, "time_series_outputs/ordination_output/score_PCoAsite_cover.csv")
write.csv(count_merge, "time_series_outputs/ordination_output/score_PCoAsite_count.csv")

write.csv(spp.score.cover, "time_series_outputs/ordination_output/spp_score_cover.csv")
write.csv(spp.score.count, "time_series_outputs/ordination_output/spp_score_count.csv")


# PLOT PCoA #####################################################################
###________________________________________________________________________________
score.noaa$year <- as.numeric(score.noaa$year)

plot_pcoa_noaa <- ggplot()+ 
  geom_point(data = score.noaa, aes(y = MDS2, x = MDS1, fill = year), shape = 21, alpha = 0.7, size = 3)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score, aes(x = MDS1, y = MDS2, label = species), max.overlaps = 20, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
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

plot_pcoa_count <- ggplot()+ 
  geom_point(data = count_merge, aes(y = MDS2, x = MDS1, fill = YEAR, shape = SITE), shape = 21, alpha = 0.7, size = 3)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score.count, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score.count, aes(x = MDS1, y = MDS2, label = species), max.overlaps = 25, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (50%)") + 
  ylab("PCoA 2 (21%)") +
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




####__________________________________________________________________
##### TEST USING NOAA SUBTIDAL SPECIES ALSO PRESENT IN THE KEEN DATA
####__________________________________________________________________


noaa_keen <- noaa_cut %>% select(year,
                                 cyclopterus_lumpus,
                                 homarus_americanus,
                                 myoxocephalus_octodecemspinosus,
                                 placopecten_magellanicus,
                                 pollachius_virens,
                                 pseudopleuronectes_americanus,
                                 tautogolabrus_adspersus,
                                 sum_cancer,
                                 sum_urophycis)

noaa_keen <- noaa_keen %>% remove_rownames %>% column_to_rownames(var="year") 


# Calculate PCoA

community <- decostand(noaa_keen, method = "hellinger") 



#PCoA
pcoa <- capscale(community~1, distance = "euclidean") # 25 et 19 %
summary(pcoa) #37%, 18%, 14%

#extraire scores sites

score.pcoa.sites <- as.data.frame(scores(pcoa, choices = c(1:3), display = "sites"))
setDT(score.pcoa.sites, keep.rownames = "year")[]
score.noaa <- score.pcoa.sites

#extraire score especes

spp.score <- data.frame(scores(pcoa, display = "species", choices = c(1:3)))
setDT(spp.score, keep.rownames = "species")[]

write.csv(score.pcoa.sites, "time_series_outputs/ordination_output/score_noaa-KEEN_pcoa_site.csv")
write.csv(spp.score, "time_series_outputs/ordination_output/score_noaa-KEEN_specie.csv")


### Run model on coordinates 

score.pcoa.sites$year <- as.numeric(score.pcoa.sites$year)
noaa_pcoa_mod <- bf( mvbind(MDS1,MDS2,MDS3) ~ s(year), family = gaussian())
get_prior(noaa_pcoa_mod, data = score.pcoa.sites)

prior_noaa_pcoa_mod <- c(prior(normal(0,1), class = "b", resp = ""),
                         prior(normal(0,1), class = "b", resp ="MDS1"),
                         prior(normal(0,1), class = "b", resp = "MDS2"),
                         prior(normal(0,1), class = "b", resp = "MDS3"),
                         prior(normal(0,1), class = "Intercept", resp = ""),
                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS1"),
                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS2"),
                         prior(student_t(3,0,2.5), class = "Intercept", resp = "MDS3"),
                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS1"),
                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS2"),
                         prior(student_t(3,0,2.5), class = "sds", resp = "MDS3"),
                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS1"),
                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS2"),
                         prior(student_t(3,0,2.5), class = "sigma", resp = "MDS3"))

## fit !!!!

fit <- brm(noaa_pcoa_mod,
           prior = prior_noaa_pcoa_mod,
           data = score.pcoa.sites,
           warmup = 2000, iter = 7000,
           cores = 3, chains = 3,
           control = list(adapt_delta = 0.99))

saveRDS(fit, file = "time_series_outputs/noaa_pcoa-KEEN_GS_3.rds")


fit
pp_check(fit, resp = "MDS1")
pp_check(fit, resp = "MDS2")
pp_check(fit, resp = "MDS3")

mcmc_plot(fit, type = "trace")
mcmc_plot(fit, type = "hist")

conditional_smooths(fit)



### plot PCoA ##################################################

score.noaa$year <- as.numeric(score.noaa$year)

plot_pcoa_noaa <- ggplot()+ 
  geom_point(data = score.noaa, aes(y = MDS2, x = MDS1, fill = year), shape = 21, alpha = 0.7, size = 3)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score, aes(x = MDS1, y = MDS2, label = species), max.overlaps = 20, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (37%)") + 
  ylab("PCoA 2 (18%)") +
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


plot_pcoa_noaa2 <- ggplot()+ 
  geom_point(data = score.noaa, aes(y = MDS3, x = MDS1, fill = year), shape = 21, alpha = 0.7, size = 3)+ #this is the points
  scale_fill_gradient(low="blue", high="red")+
  geom_segment(data = spp.score, aes(x = 0, xend = MDS1, y = 0, yend = MDS3), alpha = 0.6, #fl?ches
               arrow = arrow(length = unit(0.015, "npc"), type = "open"),
               lwd = 1)+
  geom_hline(yintercept = 0, lty = 2) + #lignes ? 0,0
  geom_vline(xintercept = 0, lty = 2) +
  geom_text_repel(data = spp.score, aes(x = MDS1, y = MDS3, label = species), max.overlaps = 20, force = 1.1, point.padding = 0.1, hjust = T, size = 4.8)+   ## text esp?ces dans NMDS
  xlab("PCoA 1 (37%)") + 
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
plot_pcoa_noaa2

ggsave("time_series_outputs/figures/supplement/pcoa_noaa-KEEN-test.tiff",
       plot = plot_pcoa_noaa, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)


ggsave("time_series_outputs/figures/supplement/pcoa_noaa-KEEN-test2.tiff",
       plot = plot_pcoa_noaa2, device = "tiff", scale = 2, 
       width = 1920, height = 991, units = "px", dpi = 300, 
       bg = "white",
       compression = "lzw"
)