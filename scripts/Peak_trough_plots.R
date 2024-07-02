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

### ###### ###### ###### ###### ###### ###### ###
######## Plots for BES and evolution talk ### ###
### ###### ###### ###### ###### ###### ###### ###
library(ggrepel)
cbPalette <- c("#F0E442","#D55E00" , "#0072B2","#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black", "#999999")

p.BES <- ggplot(raw.pheno[!is.na(raw.pheno$CLUSTER),]) +
  geom_point(aes(Julian.date, asin(sqrt(draft.estimate.prop.pigment)), col = cluster_season), alpha = 0.5) +
  geom_segment(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),],
               aes(x = time.plot, xend = time.plot, y = mean-se, yend = mean+se, col = time), size = 2, show.legend = F) +
  geom_point(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),], 
             aes(x = time.plot, y = mean, fill = time), col = "black", shape = 23, size = 5, show.legend = F, stroke = 1.5) +
  geom_text_repel(data = peak.trough[peak.trough$CLUSTER!="all data"&!is.na(peak.trough$CLUSTER),], 
             aes(x = time.plot, y = mean, label = CLUSTER), size = 8, col = "black", min.segment.length = 0, show.legend = F, nudge_x = 20) +
  scale_color_manual(values = c("grey80", "darkred", "deepskyblue"), labels = c("off-season", "peak", "trough")) +
  scale_fill_manual(values = c("darkred", "deepskyblue"), labels = c("off-season", "peak", "trough")) +
  ylim(c(asin.limits)) +
  ylab("Wing melanisation (arc-sine & square-root)") +
  xlab("Day of Year") +
  theme_bw() +
  labs(color = "Season") +
  theme(text = element_text(size = 20)) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(override.aes = list(size=6, alpha = 1)))


p.r.t <- p.BES + 
  facet_wrap(~CLUSTER, ncol = 1, strip.position = "right") +
  theme(panel.background = element_rect(fill='white'), #transparent panel bg
      plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
      #panel.grid.major = element_blank(), #remove major gridlines
      #panel.grid.minor = element_blank(), #remove minor gridlines
      legend.box.background = element_rect(fill='transparent'),
      legend.background = element_rect(fill='transparent', color = '#E5D9BA'),
                                       legend.key = element_blank()
      ) #transparent legend panel)
ggsave(filename = "plots/pheno_plot_transparent.png", plot = p.r.t, width = 12, height = 9)

p.r.t <- p.BES + theme(panel.background = element_rect(fill='white'), #transparent panel bg
                     plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
                     #panel.grid.major = element_blank(), #remove major gridlines
                     #panel.grid.minor = element_blank(), #remove minor gridlines
                     legend.box.background = element_rect(fill='transparent'),
                     legend.background = element_rect(fill='transparent', color = '#E5D9BA'),
                     legend.key = element_blank()) #transparent legend panel)
ggsave(filename = "plots/pheno_plot_transparent_single plot.png", plot = p.r.t, width = 8, height = 7)

# Map Plot
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
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(text = element_text(size = 20)) +
  theme(plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
        #panel.grid.major = element_blank(), #remove major gridlines
        #panel.grid.minor = element_blank(), #remove minor gridlines
        legend.box.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent', color = '#E5D9BA'))
##theme(text = element_text(size = 20), axis.text = element_blank(), axis.title = element_blank(),
##      axis.ticks = element_blank(),
##      plot.background = element_rect(fill=rgb(1,1,1, alpha = 0.5)),
##      plot.margin = unit(c(0,0,0,0), units = "mm")) 
map.plot  
ggsave(filename = "plots/map_plot_transparent_single plot.png", plot = map.plot  , width = 8, height = 8)

tree.plot.final.t <- tree.plot.final + theme(panel.background = element_rect(fill = 'transparent'),
                                                panel.ontop = TRUE) +
  theme(text = element_text(size = 20)) +
  theme(plot.background = element_rect(fill='#E5D9BA', color='transparent'))

ggsave(filename = "plots/tree_plot_transparent_single plot.png", plot = tree.plot.final.t  , width = 5, height = 5)


## Simulation examples
# Plot simulated data and true data
# Trait values with row names corrisponding to tips
cbPalette <- c("#F0E442","#D55E00" , "#0072B2","#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black", "#999999")

var.mean <- cbind.data.frame(peak = peak.trough$mean[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)], 
                             trough = peak.trough$mean[peak.trough$time=="trough"&!is.na(peak.trough$CLUSTER)])
row.names(var.mean) <- peak.trough$CLUSTER[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)]

# Errors on trait values
var.se <- cbind.data.frame(peak = peak.trough$se[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)], 
                           trough = peak.trough$se[peak.trough$time=="trough"&!is.na(peak.trough$CLUSTER)])
row.names(var.se) <- peak.trough$CLUSTER[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)]
# Convert standar error into squared standard errors
var.sq.se <- var.se^2

#### TREES ####
no.ce.matrix <- matrix(c(1,NA,NA,1), nrow = 2)
# create tree with different regemes for Pacfic and Atlantic
# Separtate Pacific and Atlantic into different regemes
tree.phylo.r1 <- paintSubTree(tree.phylo, node = 9, state = "All")
tree.phylo.r2.P <- paintSubTree(tree.phylo, node=10, state="Pacific", anc.state="group_1",stem=TRUE)
# Separate most northern lineage 
tree.phylo.r2.HL <- paintSubTree(tree.phylo, node=6, state="High_Lat", anc.state="group_1",stem=TRUE)
# Separate both Pacific and most northern
tree.phylo.r3.HL.P <- paintSubTree(tree.phylo.r2.HL, node=10, state="Pacific", anc.state="group_1",stem=TRUE)

# Test for different selective regimes
# Model brownian motion
nsims=1000
asin.limits <- c(asin(sqrt(c(0,1))))
modelOUr1.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1)
modelOUr1.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1,
                        param = list(sigma="constraint", alpha="constraint"))
# Model OU with and without co-evolution


modelOUr1.no.ce.sims <- mvSIM(tree = tree.phylo, nsim = nsims, model = c("OU1"), param = list(theta = modelOUr1.no.ce$theta,
                                                                                              alpha = modelOUr1.no.ce$alpha,
                                                                                              sigma = modelOUr1.no.ce$sigma))

# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
lapply(modelOUr1.no.ce.sims, function(x) points(x, col = as.factor(row.names(x))))
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

modelOUr1.no.ce.sims <- as.data.frame(do.call("rbind",modelOUr1.no.ce.sims))
modelOUr1.no.ce.sims$cluster <- stringr::str_split_i(rownames(modelOUr1.no.ce.sims),"\\.",1)
head(modelOUr1.no.ce.sims)
modelOUr1.no.ce.sims$oobs <- modelOUr1.no.ce.sims$V1<asin.limits[1]|modelOUr1.no.ce.sims$V1>asin.limits[2]|modelOUr1.no.ce.sims$V2<asin.limits[1]|modelOUr1.no.ce.sims$V2>asin.limits[2]

modelOUr1.no.ce.sims.plot <- ggplot(modelOUr1.no.ce.sims[!modelOUr1.no.ce.sims$oobs,]) +
  geom_point(aes(V1,V2,col = cluster), alpha = 0.5, show.legend = F) +
  geom_point(data = var.mean, aes(peak, trough, fill = LETTERS[1:8]), col = "black", size = 5, shape = 21, stroke = 2) +
  scale_color_manual(values = cbPalette) + scale_fill_manual(values = cbPalette) +
  labs(x = "peak", y = "off-peak") +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(text = element_text(size = 18)) +
  theme(plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
        #panel.grid.major = element_blank(), #remove major gridlines
        #panel.grid.minor = element_blank(), #remove minor gridlines
        legend.box.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent', color = '#E5D9BA'), 
        legend.key = element_blank(), legend.title = element_blank())

ggsave(filename = "plots/simulation_plot_OUr1_no_ce.png", plot = modelOUr1.no.ce.sims.plot, width = 6, height = 5)

sim.OU.r1.r3.HL.P.LRT <- readRDS("data/models/sim.OU.r1.r3.HL.P.LRT.rds")
hist(sim.OU.r1.r3.HL.P.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)

# Number of simulated LRT that were greater than observed LRT
LRT.OU.r1.r3.HL.P <- table(sim.OU.r1.r3.HL.P.LRT$LRT.lim > LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100

sim.OU.r1.r3.HL.P.LRT$exceeds.lim
t <- ggplot(sim.OU.r1.r3.HL.P.LRT[sim.OU.r1.r3.HL.P.LRT$LRT.lim>0,]) +
  geom_histogram(aes(LRT.lim,
                     fill = sim.OU.r1.r3.HL.P.LRT$exceeds.lim[sim.OU.r1.r3.HL.P.LRT$LRT.lim>0]), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, y = 100*0.8), size = 8) +
  theme_bw() +
  ylim(c(0, 100)) +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("High latitude, Pacific, and Atlantic") +
  labs(subtitle = paste0(LRT.OU.r1.r3.HL.P[2], "% greater than the observed models")) +
  theme(plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
        #panel.grid.major = element_blank(), #remove major gridlines
        #panel.grid.minor = element_blank(), #remove minor gridlines
        legend.box.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent', color = '#E5D9BA'), 
        legend.key = element_blank(), legend.title = element_blank(),
        text = element_text(size = 18))

t

modelOUr1.no.ce.sims.plot / t

ggsave(filename = "plots/simulation_plot_OUr1_no_ce_with_hist.png", plot = modelOUr1.no.ce.sims.plot / t, width = 6, height = 10)


#### Simulation without co-evolution

modelOUr3.HL.P.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P, param = list(sigma="constraint", alpha="constraint"))

model.sims.OUr3.HL.P.no.ce <- mvSIM(tree = tree.phylo.r3.HL.P, 
                                    nsim = nsims, model = c("OUM"), 
                                    param = list(theta = modelOUr3.HL.P.no.ce$theta,
                                                 alpha = modelOUr3.HL.P.no.ce$alpha,
                                                 sigma = modelOUr3.HL.P.no.ce$sigma))



# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
lapply(model.sims.OUr3.HL.P.no.ce, function(x) points(x, col = as.factor(row.names(x))))
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

model.sims.OUr3.HL.P.no.ce <- as.data.frame(do.call("rbind",model.sims.OUr3.HL.P.no.ce))
model.sims.OUr3.HL.P.no.ce$cluster <- stringr::str_split_i(rownames(model.sims.OUr3.HL.P.no.ce),"\\.",1)
model.sims.OUr3.HL.P.no.ce$oobs <- model.sims.OUr3.HL.P.no.ce$V1<asin.limits[1]|model.sims.OUr3.HL.P.no.ce$V1>asin.limits[2]|model.sims.OUr3.HL.P.no.ce$V2<asin.limits[1]|model.sims.OUr3.HL.P.no.ce$V2>asin.limits[2]

head(model.sims.OUr3.HL.P.no.ce)


model.sims.OUr3.HL.P.no.ce.plot <- ggplot(model.sims.OUr3.HL.P.no.ce[!model.sims.OUr3.HL.P.no.ce$oobs,]) +
  geom_point(aes(V1,V2,col = cluster), alpha = 0.5, show.legend = F) +
  geom_point(data = var.mean, aes(peak, trough, fill = LETTERS[1:8]), col = "black", size = 5, shape = 21, stroke = 2) +
  labs(x = "peak", y = "off-peak") +
  scale_color_manual(values = cbPalette) + scale_fill_manual(values = cbPalette) +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(text = element_text(size = 20)) +
  theme(plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
        #panel.grid.major = element_blank(), #remove major gridlines
        #panel.grid.minor = element_blank(), #remove minor gridlines
        legend.box.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent', color = '#E5D9BA'), 
        legend.key = element_blank(), legend.title = element_blank())


ggsave(filename = "plots/simulation_plot_OUr3_HL_P_no_ce.png", plot = model.sims.OUr3.HL.P.no.ce.plot, width = 6, height = 5)

sim.OU.r3.HL.P.ce.LRT <- readRDS("data/models/sim.OU.r3.HL.P.ce.LRT.rds")
hist(sim.OU.r3.HL.P.ce.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)
summary(sim.OU.r3.HL.P.ce.LRT$LRT)

LRT.OU.r3.ce <- table(sim.OU.r3.HL.P.ce.LRT$LRT.lim>LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100


k <- ggplot() +
  geom_histogram(aes(sim.OU.r3.HL.P.ce.LRT$LRT.lim, 
                     fill = as.character(modelOUr1.no.ce.sims.LRT$exceed.limits)), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio, y = (100)*0.8), size = 8) +
  theme_bw() +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  ylim(c(0, 100)) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Co-evolution with High Latitude, Pacific, and Atlantic regimes") +
  labs(subtitle = paste0(LRT.OU.r3.ce[2], "% greater than the observed models")) +
  theme(text = element_text(size = 18), title = element_text(size = 18*0.6),
        plot.background = element_rect(fill='#E5D9BA', color=NA), #transparent plot bg
        #panel.grid.major = element_blank(), #remove major gridlines
        #panel.grid.minor = element_blank(), #remove minor gridlines
        legend.box.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent', color = '#E5D9BA'), 
        legend.key = element_blank(), legend.title = element_blank())

table(sim.OU.r3.HL.P.ce.LRT$corr.lim>cov2cor(stationary(modelOUr3.HL.P))[1,2])/nsims*100
table(abs(sim.OU.r3.HL.P.ce.LRT$corr.lim)<abs(cov2cor(stationary(modelOUr3.HL.P))[1,2]))/nsims*100

kk <- ggplot() +
  geom_histogram(aes(sim.OU.r3.HL.P.ce.LRT$corr.lim, fill = as.character(modelOUr1.no.ce.sims.LRT$exceed.limits)), binwidth = 0.05, show.legend = F) + 
  geom_vline(xintercept = cov2cor(stationary(modelOUr3.HL.P))[1,2]) +
  geom_point(aes(x = cov2cor(stationary(modelOUr3.HL.P))[1,2], y = 50*0.8), size = 6, shape = 17) +
  theme_bw() +
  scale_fill_manual(values = c("orange", "orangered")) +
  ylim(c(0, 50)) +
  xlim(-1,1) +
  xlab("Standardised stationary covariance") +
  ylab("Frequency") +
  theme(text = element_text(size = 12), title = element_text(size = 12*0.6),
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        #panel.grid.major = element_blank(), #remove major gridlines
        #panel.grid.minor = element_blank(), #remove minor gridlines
        legend.box.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent', color = '#E5D9BA'), 
        legend.key = element_blank(), legend.title = element_blank())

k <- k + inset_element(kk, align_to = "plot", left = 0.5, right = 0.95, bottom = 0.4, top = 0.95, ignore_tag = F)

model.sims.OUr3.HL.P.no.ce.plot / k

ggsave(filename = "plots/simulation_plot_OUr3_no_ce_with_hist.png", plot = model.sims.OUr3.HL.P.no.ce.plot / k, width = 6, height = 10)
