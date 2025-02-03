## Get phenotype data
library(patchwork)
library(tidyverse)
library(ggrepel)

#Worldmap
world_map <- read_sf("data/world_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
e <- raster::extent(-130, -50, 0, 50)
world_map <- st_crop(world_map,e)

# read in level 3 and 4 of hydrosheds
hydrobasins_3 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev03_v1c.shp")
hydrobasins_4 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev04_v1c.shp")
hydrobasins_5 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev05_v1c.shp")

# Phenotype data
phen<-read.csv("data/MASTER_titia_phenotypes_v1.0.csv")

## Get samples that were included WGS phylogeny 
# Hydro basin 5
WGS_hydro5 <- read.table("data/SNPs/SNAPP/titia.mysnps-SNAPP-hydro_5_max_cov.txt", header = T)
# Get hydrobasin name
WGS_hydro5$hydro5 <- substr(WGS_hydro5$species, 1, 10)
#Hydrobasin4
## Get samples that were included WGS phylogeny 
WGS_hydro4 <- read.table("data/SNPs/SNAPP/titia.mysnps-SNAPP-hydro_4_max_cov.txt", header = T)
# Get hydrobasin name
WGS_hydro4$hydro4 <- substr(WGS_hydro4$species, 1, 10)

## Which phenotypes are from basins that include WGS samples
table(phen$HYBAS_ID_5%in%WGS_hydro5$hydro5)
phen$WGS_hydro_ID_5 <- phen$HYBAS_ID_5%in%WGS_hydro5$hydro5

# plot where phenotypes are and which are in the same basin as WGS samples
ggplot(phen) +
  geom_point(data = phen[!phen$WGS_hydro_ID_5,], aes(longitude, latitude), color = "grey") +
  geom_point(data = phen[phen$WGS_hydro_ID_5,], aes(longitude, latitude,col = as.factor(HYBAS_ID_5)))

# Plot all pheno data
p <- ggplot(phen[phen$WGS_hydro_ID_5,]) +
  geom_point(aes(Julian.date, draft.estimate.prop.pigment, col = as.factor(HYBAS_ID_5)), show.legend = F) +
  facet_wrap(~HYBAS_ID_5, ncol = 1)

# Plotgeographic location
q <- ggplot(phen[phen$WGS_hydro_ID_5,]) +
  geom_sf(data = world_map) +
  geom_point(aes(longitude, latitude,col = as.factor(HYBAS_ID_5)), show.legend = F) +
  facet_wrap(~HYBAS_ID_5, ncol = 1)

## Identified clusters
# Pre-lim list of identified clusters
clusters <- list()
clusters[[1]] <- c(7050000010)
clusters[[2]] <- c(7050004220, 7050005420)
clusters[[3]] <- c(7050002710)
clusters[[4]] <- c(7050047850, 7050736840, 7050047980, 7050048620, 7050047970, 7050683900)
clusters[[5]] <- c(7050668250, 7050741770, 7050682360)
clusters[[6]] <- c(7050822250, 7050823030)
clusters[[7]] <- c(7050050000, 7050050930, 7050849600)
clusters[[8]] <- c(7050050640, 7050050650)
clusters[[9]] <- c(7050054660, 7050054670)
clusters[[10]] <- c(7050043860, 7050043870, 7050045480)
clusters[[11]]  <- c(7050045900)
clusters[[12]]  <- c(7050046750)
clusters[[13]]  <- c(7050041410, 7050453940, 7050494360)
clusters[[14]]  <- c(7050609780)
clusters[[15]]  <- c(7050511950, 7050531210)

# All hydrobasins included in cluster
all_hydro <- unlist(clusters)
# Any WGS samples not used in clusters?
miss_hydro <- WGS_hydro5[which(!WGS_hydro5$hydro5%in%all_hydro),]$hydro5
# Hydrobasins that have been assigned to basin without direct WGS samples
extra_hydro <- all_hydro[!all_hydro%in%WGS_hydro5$hydro5]

p <- ggplot(phen[!is.na(phen$WGS_hydro_ID_5),]) +
  geom_sf(data = world_map) +
  geom_sf(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%all_hydro,], fill = "lightgreen", col = "black") +
  geom_sf(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%miss_hydro,], fill = "grey60", col = "black") +
  geom_sf(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%extra_hydro,], fill = "purple1", col = "black") +
  geom_sf_text(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%all_hydro,], aes(label=HYBAS_ID), size = 1) +
  geom_sf_text(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%extra_hydro,], aes(label=HYBAS_ID), size = 1) +
  geom_sf_text(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%miss_hydro,], aes(label=HYBAS_ID), size = 1)

ggsave("plots/cluster_exploration/unassigned_WGS_samples.pdf", p)

# Loop that assigns each phenotype to a sample
phen$WGS_hydro_ID_5 <- NA
for(i in 1:length(clusters)){
  phen$WGS_hydro_ID_5[phen$HYBAS_ID_5%in%clusters[[i]]] <- paste0("Clust",i)
}
# Breakdown of assignment
summary(as.factor(phen$WGS_hydro_ID_5))

# Summary of lat and long of each cluster
clust_summ <- group_by(phen, by = WGS_hydro_ID_5) %>%
  summarise(mn.lat = median(latitude), mn.lon = median(longitude))
clust_summ <- clust_summ[1:length(clusters),]

# map with all clusters
q <- ggplot(phen[!is.na(phen$WGS_hydro_ID_5),]) +
  geom_sf(data = world_map) +
  geom_sf(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%phen$HYBAS_ID_5,], fill = "grey60", col = "black") +
  geom_sf(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%all_hydro,], fill = "grey40", col = "black") +
  geom_sf_text(data = hydrobasins_5[hydrobasins_5$HYBAS_ID%in%phen$HYBAS_ID_5,], aes(label = HYBAS_ID),size=0.5) +
  geom_point(aes(longitude, latitude,col = as.factor(WGS_hydro_ID_5)), show.legend = F, size=0.1)
  #geom_text_repel(data = clust_summ, aes(mn.lon, mn.lat, label = by),min.segment.length = 0,force = 10)

p <- ggplot(phen) +
  geom_point(aes(x = Julian.date, y = draft.estimate.prop.pigment, col = as.factor(WGS_hydro_ID_5)), size=0.2) +
  facet_wrap(~WGS_hydro_ID_5)


ggsave("plots/cluster_exploration/clusters_mapped.pdf", q+p, height = 10,width = 15)


#### Assigning samples outside of basins with 