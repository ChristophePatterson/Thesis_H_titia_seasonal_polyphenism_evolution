## read in phenotype daa
source("scripts/Data_extraction.R")

library(shadowtext)

## Get peak and trough values
peak.trough$time.plot[peak.trough$time=="peak"] <- rowMeans(peak.trough.season[,2:3])
peak.trough$time.plot[peak.trough$time=="trough"] <- 330

## Plot of all data
library(ggrepel)
cbPalette <- c("#F0E442","#D55E00" , "#0072B2","#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black", "#999999")

# Calculate logit
raw.pheno$proportion.pigmented.logit <- log( (raw.pheno$proportion.pigmented-0.01) / (1-((raw.pheno$proportion.pigmented-0.01))))

# Plot1
p.pheno <- ggplot(raw.pheno[!is.na(raw.pheno$CLUSTER),]) +
  geom_point(aes(Julian.date, proportion.pigmented.logit, fill = CLUSTER, alpha = cluster_season), size = 3, shape = 21, show.legend = F) +
  geom_point(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),], 
             aes(x = time.plot, y = mean, fill = CLUSTER), col = "black", shape = 23, size = 5, stroke = 1.5) +
  geom_label_repel(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),], 
                  aes(x = time.plot, y = mean, label = CLUSTER), fill = rgb(1,1,1,0.75), size = 8, col = "black", min.segment.length = 0, show.legend = F, nudge_x = 30) +
  # scale_color_manual(values = c("grey80", "darkred", "deepskyblue"), labels = c("off-season", "peak", "trough")) +
  scale_fill_manual(values = cbPalette) +
  scale_alpha_manual(values = c(0.1,1,1)) +
  ylim(c(-5, 5)) +
  ylab("Wing melanisation (logit)") +
  xlab("Day of Year") +
  theme_bw() +
  labs(color = "Season") +
  theme(text = element_text(size = 18)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(override.aes = list(size=6, alpha = 1)))


p.pheno

library(sf)

# Geographic distibution
world_map <- read_sf("data/world_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#ice_sheet <- read_sf("data/world_data//Ice Sheets/18ka_22.1cal_Dalton_et_al_2020_QSR.shp")
e <- raster::extent(-130, -50, 0, 50)
world_map <- st_crop(world_map,e)
crs_text <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",360+mean(raw.pheno$longitude), " +lat_0=",mean(raw.pheno$latitude))
worldmap_proj <- st_transform(world_map, crs = crs_text)
plot(world_map[1])

# map of samples
all_hyrdo <- lapply(unique(raw.pheno$CLUSTER), FUN = function(x) {
  print(x)
  return(cbind(unique(raw.pheno$HYBAS_ID_5[raw.pheno$CLUSTER==x]), as.character(x)))})
all_hyrdo 
all_hyrdo <- do.call("rbind", all_hyrdo)

table(raw.pheno$CLUSTER, raw.pheno$HYBAS_ID_5)

pheno.location.summary <- raw.pheno %>%
  group_by(CLUSTER) %>%
  summarise(lat = mean(latitude, na.rm = T), lon = mean(longitude, na.rm = T))

# Find the convex hull of the points being plotted 
hull <- raw.pheno %>% group_by(CLUSTER) %>% 
  slice(chull(longitude, latitude))

map.plot <- ggplot(raw.pheno) +
  geom_sf(data = worldmap_proj, fill = "grey90", col = "black", linewidth = 0.5) +
  #geom_sf(data = ice_sheet, fill = "skyblue", alpha = 0.5) +
  geom_point(aes(longitude, latitude, col = CLUSTER), size = 2, show.legend = F) +
  #geom_polygon(data = hull[hull$CLUSTER!="NA",], aes(longitude, latitude, colour = CLUSTER,), size = 1,  fill = NA, show.legend = F) +
  #geom_text(data = pheno.location.summary, aes(lon, lat, label = CLUSTER), size = 11, col = "black") +
  #geom_shadowtext(data = pheno.location.summary, aes(lon, lat, label = CLUSTER), size = 8, bg.colour='black') +
  coord_sf(xlim = range(raw.pheno$longitude)+c(-5,+5), ylim =range(raw.pheno$latitude)+c(-2,+2), default_crs = st_crs(4326)) +
  scale_color_manual(values = cbPalette) +
  labs(x = "Longitude", y = "Latitude", col = "Tree tip") +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(text = element_text(size = 18))


map.plot

## Tree plot
##not the most elegant solution for reordering tips, but produces an accurate tree

tree.phylo$tip.label<-c("A","B","E","F","H","G","C","D")
tree.plot <- ggtree(tree.phylo, size = 1.2) 

tree.map <- (tree.plot + geom_tiplab(aes(x = x + 0.1), size = 6) + geom_tippoint(aes(fill = label), col = "black", shape = 23, size = 5, stroke = 1.5, show.legend = F, ) + scale_fill_manual(values = cbPalette, breaks = plot.tip.order) 
             +geom_treescale(width=0.5,linesize=1.2)+ scale_x_continuous(expand = c(0.1, 0.1))) +  theme(text = element_text(size = 18))

figure1 <- p.pheno / (tree.map + map.plot) / guide_area() + plot_layout(heights = c(2,1.75, 0.2), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("plots/Figure1.png", figure1, width = 7.24*1.8, height = 7.24*2)

