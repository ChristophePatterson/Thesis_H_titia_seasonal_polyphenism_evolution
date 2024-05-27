## Read in data
source("scripts/Data_extraction.R")
library(sf)
# install.packages("shadowtext")

asin.limits <- c(asin(sqrt(c(0,1))))

p <- ggplot(raw.pheno[!is.na(raw.pheno$CLUSTER),]) +
  geom_point(aes(Julian.date, asin(sqrt(draft.estimate.prop.pigment)), col = cluster_season)) +
  scale_color_manual(values = c("grey80", "darkred", "deepskyblue"), labels = c("off-season", "peak", "trough")) +
  ylim(c(asin.limits)) +
  facet_wrap(~CLUSTER, ncol = 1, strip.position = "right") +
  ylab("Wing melanisation (arc-sine & square-root)") +
  xlab("Day of Year") +
  theme_bw() +
  labs(color = "Season") +
  theme(text = element_text(size = 20)) +
  theme(legend.position="bottom")

peak.trough$time.plot[peak.trough$time=="peak"] <- c(NA,rowMeans(peak.trough.season[,2:3]))
peak.trough$time.plot[peak.trough$time=="trough"] <- 330

p.r <- ggplot(raw.pheno[!is.na(raw.pheno$CLUSTER),]) +
  geom_point(aes(Julian.date, asin(sqrt(draft.estimate.prop.pigment)), col = cluster_season), alpha = 0.5) +
  geom_segment(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),],
               aes(x = time.plot, xend = time.plot, y = mean-se, yend = mean+se, col = time), size = 2, show.legend = F) +
  geom_point(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),], 
             aes(x = time.plot, y = mean, fill = time), col = "black", shape = 23, size = 5, show.legend = F, stroke = 1.5) +
  scale_color_manual(values = c("grey80", "darkred", "deepskyblue"), labels = c("off-season", "peak", "trough")) +
  scale_fill_manual(values = c("darkred", "deepskyblue"), labels = c("off-season", "peak", "trough")) +
  ylim(c(asin.limits)) +
  facet_wrap(~CLUSTER, ncol = 1, strip.position = "right") +
  ylab("Wing melanisation (arc-sine & square-root)") +
  xlab("Day of Year") +
  theme_bw() +
  labs(color = "Season") +
  theme(text = element_text(size = 20)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(override.aes = list(size=6, alpha = 1)))

q <- ggplot(peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),]) +
  geom_point(aes(x = time, y = mean, col = time), size = 3,show.legend = F) +
  geom_segment(aes(x = time, xend = time, y = mean-se, yend = mean+se, col = time), size = 2, show.legend = F) +
  scale_color_manual(values = c("darkred", "deepskyblue")) +
  ylim(asin.limits) +
  facet_wrap(~CLUSTER, ncol = 1, strip.position = "right") +
  ylab("Wing melanisation (arc-sine & square-root)") +
  xlab("Season") +
  labs(color = "Season") +
  theme_bw() +
  theme(text = element_text(size = 20))

# How old are the splits between titia
(((max(tree.plot$data$x)/tree.plot$data$x)*3.5)-3.5)*1000000
tree.plot + geom_label(aes(label = round((x-max(x))*1000000, 2)))

tree.plot$data$x


tree.plot.final <- tree.plot + theme(text = element_text(size = 18)) + xlim(min(tree.plot$data$x), max(tree.plot$data$x)+0.1)
tree.plot.final
# Geographic distibution
world_map <- read_sf("data/world_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
ice_sheet <- read_sf("data/world_data//Ice Sheets/18ka_22.1cal_Dalton_et_al_2020_QSR.shp")
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

cbPalette <- c("#F0E442","#D55E00" , "#0072B2","#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black", "#999999")

pheno.location.summary <- raw.pheno %>%
  group_by(CLUSTER) %>%
  summarise(lat = mean(latitude, na.rm = T), lon = mean(longitude, na.rm = T))

# Find the convex hull of the points being plotted 
hull <- raw.pheno %>% group_by(CLUSTER) %>% 
  slice(chull(longitude, latitude))

ggplot(raw.pheno) +
  geom_point(aes(longitude, latitude, colour = CLUSTER)) +
  geom_polygon(data = hull[hull$CLUSTER!="NA",], alpha = 0.2, 
                    aes(longitude, latitude, fill = CLUSTER, colour = CLUSTER))


library(shadowtext)
map.plot <- ggplot(raw.pheno) +
  geom_sf(data = worldmap_proj, fill = "grey90", col = "black", linewidth = 0.5) +
  #geom_sf(data = ice_sheet, fill = "skyblue", alpha = 0.5) +
  geom_point(aes(longitude, latitude, col = CLUSTER), show.legend = F) +
  geom_polygon(data = hull[hull$CLUSTER!="NA",], aes(longitude, latitude, colour = CLUSTER,), size = 1,  fill = NA, show.legend = F) +
  #geom_text(data = pheno.location.summary, aes(lon, lat, label = CLUSTER), size = 11, col = "black") +
  geom_shadowtext(data = pheno.location.summary, aes(lon, lat, label = CLUSTER), size = 8, bg.colour='black') +
  coord_sf(xlim = range(raw.pheno$longitude)+c(-5,+5), ylim =range(raw.pheno$latitude)+c(-2,+2), default_crs = st_crs(4326)) +
  scale_color_manual(values = cbPalette) +
  labs(x = "Longitude", y = "Latitude", col = "Tree tip") +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(text = element_text(size = 20))
  ##theme(text = element_text(size = 20), axis.text = element_blank(), axis.title = element_blank(),
  ##      axis.ticks = element_blank(),
  ##      plot.background = element_rect(fill=rgb(1,1,1, alpha = 0.5)),
  ##      plot.margin = unit(c(0,0,0,0), units = "mm")) 
map.plot  
# Save to plot
pqt <- tree.plot.final + 
  p + q + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",  theme = theme(plot.title = element_text(size = 50))) +
  plot_layout(widths = c(2,2,1))
pqt

pq <- tree.plot.final + 
  p.r + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",  theme = theme(plot.title = element_text(size = 50))) +
  plot_layout(widths = c(1,1))
pq


pqtm  <- map.plot + tree.plot.final + p.r + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.title = element_text(size = 50))) +
  plot_layout(nrow = 1, width = c(3,1.8,3))

pqtm  <- (tree.plot.final + inset_element(map.plot, left = 0, right = 0.7, top = 1, bottom = 0.58,ignore_tag = T)) + p.r + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.title = element_text(size = 50))) +
  plot_layout(nrow = 1, width = c(1.5, 1))

pqtm_test  <- (tree.plot.final / map.plot ) | p.r

pqtm_test <- pqtm_test +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")", theme = theme(plot.title = element_text(size = 50))) 
#plot_layout(height = c(2,1), width = c(1,1))


ggsave(filename = "plots/peak_trough.png", plot = pqt, width = 15, height = 15)
ggsave(filename = "plots/peak_trough_reduced.png", plot = pq, width = 15, height = 15)
ggsave(filename = "plots/peak_trough_map.png", plot = pqtm, width = 15, height = 15)

ggsave(filename = "plots/peak_trough_map_test.png", plot = pqtm_test, width = 15, height = 15)


tree.map <- (tree.plot.final + geom_tippoint(size = 3, aes(col = label), show.legend = F) + scale_color_manual(values = cbPalette)) + map.plot

ggsave(filename = "plots/tree_map_pheno.png", plot = tree.map, width = 10, height = 7.5)
ggsave(filename = "plots/pheno_map.png", plot = map.plot, width = 12, height = 7.5)


### Plot of hybrid region
# Hybrid region
world_map_hyb <- read_sf("data/world_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
load("data/world_data/hydroRIVERS/HydroRIVERS_v10_na_sub_hybrid.rda")
e.hyb <- raster::extent(-95.965576,-93.861694,15.998295,17.395200)
world_map_hyb <- st_crop(world_map_hyb,e.hyb)
plot(world_map_hyb[1])
crs
raw.pheno.hyb <- raw.pheno
raw.pheno.hyb$hybrid.region <- raw.pheno$latitude<17.395200&raw.pheno$latitude>15.998295&raw.pheno$longitude<(-93.861694)&raw.pheno$longitude>(-95.965576)

ggplot() +
  geom_sf(data = world_map_hyb) +
  geom_sf(data = hydrorivers.sub) +
  geom_point(data = raw.pheno.hyb, aes(longitude, latitude, col = CLUSTER))

ggplot(raw.pheno.hyb[raw.pheno.hyb$CLUSTER=="D",]) +
  geom_point(aes(Julian.date, asin(sqrt(draft.estimate.prop.pigment)),col = hybrid.region), )
