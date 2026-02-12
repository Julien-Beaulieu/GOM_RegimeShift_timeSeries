### This script is to make binary connectance webs connecting the subtidal 
### species to the intertidal species. The interaction data are pulled from GLOBI and 
### prepared in the script dataPrep_NOAA-intertidal_web.R.
### Contributor: Julien Beaulieu

dir_int <- read.csv("data/direct_interaction_int-sub.csv") #direct interactions from subtidal to the intertidal
indir_dir_int <- read.csv("data/indirect_direct_interaction_int-sub.csv") # includes subtidal species interacting with 
                                                                          # subtidal ps. feeding on intertiadal sp.
sp_origin <- read.csv("data/species_origin_network.csv")

<<<<<<< HEAD
=======
cross_corr <- read.csv("time_series_outputs/cross_correlation_values_inter-sub.csv")

>>>>>>> 00b592a4b0d87326e006fb045c547976cde7a4f8
### Load packages

library(igraph)
library(dplyr)
library(tidyverse)
library(ggnetwork)
library(ggplot2)
library(NetIndices)
library(rglobi)

source("Scripts/function_genusSp_2_GSp.R")

##### Direct interactions ###############################################
#####____________________________________________________________________

### remove subtidal species with trophic level = 1

sp_to_rm = c(
  "Bathypolypus arcticus",
  "Crangon septemspinosa",
  "Crustacea shrimp",
  "Dichelopandalus leptocerus",
  "Dipturus laevis",
  "Geryon quinquedens",
  "Hemitripterus americanus",
  "Illex illecebrosus",
  "Lithodes maja",
  "Lophius americanus",
  "Lumpenus lumpretaeformis",
  "Lumpenus maculatus",
  "Lycenchelys verrillii",
  "Macrozoarces americanus",
  "Maurolicus weitzmani",
  "Melanostigma atlanticum",
  "Myxine glutinosa",
  "Pandalus montagui",
  "Pandalus propinquus",
  "Paralepidae",
  "Petromyzon marinus",
  "Placopecten magellanicus",
  "Pontophilus norvegicus",
  "Scomberesox saurus",
  "Spirontocaris liljeborgii",
  "Stoloteuthis leucoptera",
  "Argentina striata",
  "Urophycis chesteri",
  "Chanceon quinquedens"
)


dir_int <- dir_int %>% filter(!(source %in% sp_to_rm)) %>% 
  filter(!(target %in% sp_to_rm))


sp_origin <- sp_origin %>% filter(!(specie %in% sp_to_rm))

<<<<<<< HEAD
=======
write.csv(dir_int, "data/species_origin_directInter.csv",
          row.names = FALSE)


### make a list of correlation to weight the edges


### ID species with negative cross-correlation
# 
# edge_attr <- cross_corr %>% unite("ID", c(sub_sp, inter_sp, site), remove = FALSE) %>% # take the max lag
#   group_by(ID) %>% mutate(Correlation = max(correlation)) %>% ungroup() %>%
#   dplyr::select(-c(ID,lag, correlation)) %>% 
#   unique() %>% unite("inter_sp_site", c(inter_sp, site), sep = "-") %>% 
#   filter(Correlation < 0)  %>%
#   separate(inter_sp_site, c("inter_sp","site"), sep = "-") %>% 
#   dplyr::select(-c(Correlation,site)) %>% 
#   unique() %>% unite("inter", c(sub_sp, inter_sp), sep = "-", remove = FALSE)



>>>>>>> 00b592a4b0d87326e006fb045c547976cde7a4f8
### put the data in interaction matrices format ########

sp_list <- unique(sp_origin$specie) 

mat <- as.data.frame(matrix(nrow = length(sp_list), ncol = length(sp_list)+1)) %>% 
  mutate(V1 = sp_list) %>% column_to_rownames(var="V1")
colnames(mat) <- sp_list

# fill the matrix with 1 if there is an interaction and 0 if not
for (i in c(1:nrow(mat))) {
  for (j in c(1:ncol(mat))) {
    temp <- dir_int %>% filter(source == row.names(mat)[i] & target == colnames(mat)[j])
    if(nrow(temp) > 0){mat[i,j] = 1}
  }
}

mat[is.na(mat)] <- 0
mat_t <- as.matrix(mat)
<<<<<<< HEAD
#as.matrix(t(mat)) # convert to matrix because can't take df
=======
>>>>>>> 00b592a4b0d87326e006fb045c547976cde7a4f8


### make a lit of subtidal species that have a trophic level = 0

sp_origin %>% mutate(TrophInd(Flow = mat_t)) %>% 
  filter(origin == "subidal") %>% filter(TL < 1.000001) %>% 
  dplyr::select(specie) %>% dput()


### make an interaction graph ################################################
####__________________________________________________________________________

<<<<<<< HEAD
=======


>>>>>>> 00b592a4b0d87326e006fb045c547976cde7a4f8
### prep species origin for layout

layout_sp <- sp_origin %>% filter(specie %in% sp_list)

# make layout with trophic level from network
layout <- as.data.frame(sp_list) %>% 
  mutate(y = TrophInd(Flow = mat_t)[,1]) %>%
  mutate(x = runif(length(sp_list), min = -1, max = 1)) %>%
  dplyr::select(c(x,y)) %>% 
  as.matrix()

sub_sp <- sp_origin %>% filter(origin == "subidal") #make a list of subtidal sp.

net <- graph_from_adjacency_matrix(mat_t)
ggnet <- ggnetwork(net, layout = layout) %>% 
  mutate(zone = ifelse(name %in% sub_sp$specie, "Subtidal", "Intertidal")) %>%
  genusSp2gSp(column = "name")
  

plot <- ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend))+
<<<<<<< HEAD
  geom_edges(aes(color = zone), alpha = 0.2) +
  geom_nodetext_repel(aes(label = name, color = zone),
                fontface = "italic", size = 3) +
  theme_blank()

plot

=======
  geom_edges(aes(color = zone), alpha = 0.1) +
  geom_nodetext_repel(aes(label = name, color = zone),
                fontface = "italic", size = 2.5,
                max.overlaps = 1000) +
  scale_color_manual(values = c("darkblue","darkred"))+
  theme_blank()+
  theme(legend.position = "none")

plot


plot2 <- ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color = zone), alpha = 0.1) +
  geom_nodes(aes(color = zone), alpha = 0.7, size = 4) +
  scale_color_manual(values = c("darkblue","darkred"))+
  theme_blank()+
  theme(legend.position = "none")

plot2

>>>>>>> 00b592a4b0d87326e006fb045c547976cde7a4f8
ggsave(filename = "time_series_outputs/figures/network_sub_inter.tiff", 
       plot = plot, device = "tiff", scale = 2,
       width = 1920, height = 991, bg = "white", 
       units = "px", dpi = 300,
       compression = "lzw")



