#Hydrobasins test map
#install.packages("ggplot2")
#install.packages("elevatr")
library(sf)
library(ggplot2)
library(elevatr)
#install.packages("elevatr")
library(terra)
# options("sp_evolution_status" = 2)
worldmap <- read_sf("World/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#Cropping the extent of worldmap to reduce plotting time
e <- c(xmin = -130, xmax = -60, ymin = 0, ymax = 45)
worldmap <- st_crop(worldmap, e)

#hydrobasins <- sf::st_read("World/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev03_v1c.shp")
#hydrobasins.sub <- hydrobasins[hydrobasins$ENDO==0,]
#hydrobasins_geo <- st_geometry(hydrobasins.sub)
#plot(hydrobasins_geo)
#hydrorivers <- sf::st_read("World/hydrosheds/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp")
#hydrorivers.sub <- hydrorivers
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ORD_CLAS<=1,]
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ORD_FLOW<=5,]
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ENDORHEIC==0,]
#
# names(hydrorivers)

#mex.hyb.e <- extent(cbind(hybrid.sites$Long,hybrid.sites$Lat))
#mex.hyb.e[1] <- mex.hyb.e[1]-0.5
#mex.hyb.e[2] <- mex.hyb.e[2]+0.5
#mex.hyb.e[3] <- mex.hyb.e[3]-0.5
#mex.hyb.e[4] <- mex.hyb.e[4]+0.5
#
#mex.hyb.e <- st_as_sf(as(mex.hyb.e,"SpatialPolygons"))
#st_crs(mex.hyb.e) <- st_crs(worldmap)
#
#hydrorivers.sub <- hydrorivers
#st_crs(hydrorivers.sub) <- st_crs(worldmap)
#
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ORD_CLAS<=6,]
#hydrorivers.sub<- hydrorivers.sub[hydrorivers.sub$ORD_FLOW<=6,]
#hydrorivers.sub.mex <- st_intersection(hydrorivers.sub, mex.hyb.e)
#hydrorivers.sub.mex
#plot(hydrorivers.sub[13])
#plot(hydrorivers.sub.mex[13])

#plot(hydrorivers.sub)
#
#save(hydrorivers.sub.mex, file = "World/hydrosheds/hydroRIVERS/HydroRIVERS_v10_na_sub_hybrid.rda")
#save(hydrorivers.sub, file = "World/hydrosheds/hydroRIVERS/HydroRIVERS_v10_na_sub1.rda")
#save(hydrobasins.sub, file = "World/hydrosheds/hydroBASIN/hybas_lake_na_lev03_v1c_sub1_rda")
load("World/hydrosheds/hydroRIVERS/HydroRIVERS_v10_na_sub1.rda")
load("World/hydrosheds/hydroBASIN/hybas_lake_na_lev03_v1c_sub1_rda")
hydrorivers_geo <- st_geometry(hydrorivers.sub)
#hydrobasins.sub <- hydrobasins[hydrobasins$AREA_SQKM>20,]
hydrobasins_geo <- st_geometry(hydrobasins.sub)

#Projection
#worldmap <- rgdal::readOGR("World/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")

#crs(worldmap)
projcrs <- "+proj=longlat +datum=WGS84 +no_defs"

#Setting extent of points that cover Mexico and CR site
World <- st_as_sf(world.e,coords = c("Long","Lat"), crs = projcrs)
e <- st_bbox(World)
Mex <- st_as_sf(mex.e,coords = c("Long","Lat"), crs = projcrs)
e.Mex <- st_bbox(Mex)
Mex.z <- st_as_sf(mex.e.zoom,coords = c("Long","Lat"), crs = projcrs)
CR <- st_as_sf(CR.e,coords = c("Long","Lat"), crs = projcrs)
e.CR <- st_bbox(CR)

#Downloading elevational data
elevation <- rast(get_elev_raster(World, z = 4, prj = projcrs))
elevation.Mex <- rast(get_elev_raster(Mex, z = 6, prj = projcrs))
elevation.Mex.z <- rast(get_elev_raster(Mex.z, z = 8, prj = projcrs))
elevation.CR <- rast(get_elev_raster(CR, z = 7, prj = projcrs))
#Setting water height to 0
elevation[elevation<=0] <- 0
elevation.CR[elevation.CR<0] <- 0
elevation.Mex[elevation.Mex<0] <- 0
elevation.Mex.z[elevation.Mex.z<0] <- 0

#Calculating aspect for America
slope.raster <- terrain(elevation*10, v='slope',unit = "radians")
aspect.raster <- terrain(elevation*10, v='aspect',unit = "radians")
class(aspect.raster)
hill.raster <- shade(slope = slope.raster, aspect = aspect.raster, angle = 45, direction =  90)
#plot(hill.raster)
hill.m <- as.points(hill.raster)
hill.df <-  terra::as.data.frame(x = hill.m,geom = "XY")
head(hill.df)
colnames(hill.df) <- c("hill", "lon", "lat")

#Calculating aspect for Mexico
slope.raster.Mex <- terrain(elevation.Mex*10, v='slope',unit = "radians")
aspect.raster.Mex <- terrain(elevation.Mex*10, v='aspect',unit = "radians")
hill.raster.Mex <- shade(slope = slope.raster.Mex, aspect = aspect.raster.Mex)
hill.m.Mex  <- as.points(hill.raster.Mex)
hill.df.Mex <-  terra::as.data.frame(x = hill.m.Mex,geom = "XY")
colnames(hill.df.Mex) <- c("hill", "lon", "lat")

#Calculating aspect for Mexico zoom
slope.raster.Mex.z <- terrain(elevation.Mex.z*10, v='slope',unit = "radians")
aspect.raster.Mex.z <- terrain(elevation.Mex.z*10, v='aspect',unit = "radians")
hill.raster.Mex.z <- shade(aspect = aspect.raster.Mex.z, slope = slope.raster.Mex.z, direction = 0, angle = 45)
#plot(hill.raster.Mex.z)
hill.m.Mex.z  <- as.points(hill.raster.Mex.z)
hill.df.Mex.z <-  terra::as.data.frame(x = hill.m.Mex.z,geom = "XY")
colnames(hill.df.Mex.z) <- c("hill", "lon", "lat")

#Calculating aspect for CR
slope.raster.CR <- terrain(elevation.CR*10, v='slope',unit = "radians")
aspect.raster.CR <- terrain(elevation.CR*10, v='aspect',unit = "radians")
hill.raster.CR <- shade(slope.raster.CR, aspect.raster.CR, 45, 90)
#plot(hill.raster.CR)
hill.m.CR <- as.points(hill.raster.CR)
hill.df.CR <-  terra::as.data.frame(hill.m.CR,geom = "XY")
colnames(hill.df.CR) <- c("hill","lon", "lat")

#head(hill.df)
#pdf("test_hydrobasins.pdf")
#ggplot() +
#  geom_raster(data = hill.df.Mex, aes(lon, lat, fill = hill), alpha = 1) +
#  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4"))
#  geom_sf(data = hydrobasins_geo, col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
#  geom_sf(data = hydrorivers_geo, col = "black", lineend = "round") +
#  geom_polygon(data = worldmap, aes(long, lat, group = group), col = "black", fill = NA) +
#  theme(legend.position="none") +
#  scale_fill_gradientn(colours=c("black","white"))
#  coord_sf(xlim = c(-97,-92), ylim = c(14,19))
#dev.off()
